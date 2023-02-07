#include "arguments.H"
#include "vector.H"
#include "block.H"
#include <chrono>
#include <fftw3.h>
#include "decomposer.H"

using namespace Foam;
using namespace briscola;

// Check the pressure field including ghost cells

bool check(scalarBlock& p, scalarBlock& f, labelVector BC)
{
    labelVector N(p.shape());

    // scalar f;

    for (int i = 0; i < N.x(); i++)
    {
        for (int j = 0; j < N.y(); j++)
        {
            for (int k = 0; k < N.z(); k++)
            {
                // Ghost cell values for i,j,k = 0, N-1
                scalar px0 = 0, pxN = 0, py0 = 0, pyN = 0, pz0 = 0, pzN = 0;

                switch (BC.x())
                {
                    case 1:
                        px0 = - p(i,j,k);
                        pxN = - p(i,j,k);
                        break;

                    case 2:
                        px0 = p(i,j,k);
                        pxN = p(i,j,k);
                        break;

                    case 3:
                        px0 = - p(i,j,k);
                        pxN = p(i,j,k);
                        break;

                    case 4:
                        px0 = p(i,j,k);
                        pxN = -p(i,j,k);
                        break;

                    case 5:
                        px0 = p(N.x()-1,j,k);
                        pxN = p(0,j,k);
                        break;
                }

                switch (BC.y())
                {
                    case 1:
                        py0 = - p(i,j,k);
                        pyN = - p(i,j,k);
                        break;

                    case 2:
                        py0 = p(i,j,k);
                        pyN = p(i,j,k);
                        break;

                    case 3:
                        py0 = - p(i,j,k);
                        pyN = p(i,j,k);
                        break;

                    case 4:
                        py0 = p(i,j,k);
                        pyN = -p(i,j,k);
                        break;

                    case 5:
                        py0 = p(i,N.y()-1,k);
                        pyN = p(i,0,k);
                        break;
                }

                switch (BC.z())
                {
                    case 1:
                        pz0 = - p(i,j,k);
                        pzN = - p(i,j,k);
                        break;

                    case 2:
                        pz0 = p(i,j,k);
                        pzN = p(i,j,k);
                        break;

                    case 3:
                        pz0 = - p(i,j,k);
                        pzN = p(i,j,k);
                        break;

                    case 4:
                        pz0 = p(i,j,k);
                        pzN = -p(i,j,k);
                        break;

                    case 5:
                        pz0 = p(i,j,N.z()-1);
                        pzN = p(i,j,0);
                        break;
                }

                scalar residual =
                  ((i == 0) ?          px0 : p(i-1,j,k))
                + ((i == N.x() - 1) ?  pxN : p(i+1,j,k))
                + ((j == 0) ?          py0 : p(i,j-1,k))
                + ((j == N.y() - 1) ?  pyN : p(i,j+1,k))
                + ((k == 0) ?          pz0 : p(i,j,k-1))
                + ((k == N.z() - 1) ?  pzN : p(i,j,k+1))
                - 6.0 * p(i,j,k) - f(i,j,k);

                if(abs(residual) > 1e-10)
                {
                    FatalError
                        << "Test failed. Residual =  " << residual
                        << " at index " << labelVector(i,j,k) << endl
                        << abort(FatalError);
                }
            }
        }
    }

    return true;
}

// Tridiagonal sovler with upper and lower diagonals filled with 1's

void tridiagonalSolve(label n, double D[], double f[], double p[])
{
    List<scalar> Dh(n);
    List<scalar> fh(n);

    // Copy D and f arrays so as not to overwrite them
    for(int i=0; i<n; i++ )
    {
        Dh[i] = D[i];
        fh[i] = f[i];
    }

    // Forward substitution
    for (int i = 1; i < n; i++ )
    {
        scalar m = 1.0 / Dh[i-1];
        Dh[i] -= m;
        fh[i] -= m * fh[i-1];
    }

    // Backward substitution
    p[n-1] = fh[n-1] / Dh[n-1];
    for (int i = n-2; i >= 0; i-- )
    {
        p[i] = (fh[i] - p[i+1]) / Dh[i];
    }
}

// Cyclic variant of tridiagonalSolver()

void cyclicTridiagonalSolve(label n, double D[], double f[], double p[])
{
    List<scalar> u(n, Zero);
    tridiagonalSolve((n-1), &D[1], &f[1], &u[1]);

    List<scalar> v(n, Zero);
    List<scalar> vf(n, Zero);

    vf[1] = -1;
    vf[n-1] = -1;

    tridiagonalSolve((n-1), &D[1], &vf[1], &v[1]);

    p[0] = (f[0]-u[n-1]-u[1]) / (D[0]+v[n-1]+v[1]);

    for (int k=1; k<n; k++ )
    {
        p[k] = u[k] + p[0]*v[k];
    }
}

int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    using std::rand;
    using std::srand;

    using Foam::pow;
    using Foam::sin;

    arguments::addBoolOption("parallel", "run in parallel");

    arguments::validArgs.append("mesh size");
    arguments::validArgs.append("initial decomposition");
    arguments::validArgs.append("boundary conditions");

    arguments args(argc, argv);

    const labelVector N(args.argRead<labelVector>(1));
    const labelVector I(args.argRead<labelVector>(2));
    const labelVector BC(args.argRead<labelVector>(3));

    // Some checks

    if (!Pstream::parRun())
    {
        FatalError << "Must be run in parallel" << endl;
        FatalError.exit();
    }

    int rank = Pstream::myProcNo();

    // Initialize decomposer

    autoPtr<Decomposer> decomp;
    decomp = new Decomposer(N, I);

    List<labelVector> Ni = decomp->Ni_;
    List<labelVector> si = decomp->si_;

    // FFT variables

    int fft_rank, howmany, stride, sep;

    fftw_plan plan_fwd_x;
    fftw_plan plan_bwd_x;

    fftw_plan plan_fwd_y;
    fftw_plan plan_bwd_y;

    scalar normalization;
    if( BC.x() != 5 && BC.y() != 5)
    {
        normalization = 2.0 * N.x() * 2.0 * N.y();
    }
    else if ( BC.x() == 5 && BC.y() == 5 )
    {
        normalization = N.x() * N.y();
    }
    else
    {
        normalization = 2.0 * N.x() * N.y();
    }

    fftw_r2r_kind transforms[] =
    {
        // forward transform types
        FFTW_RODFT10, FFTW_REDFT10, FFTW_RODFT11, FFTW_REDFT11, FFTW_R2HC,
        // inverse transform types
        FFTW_RODFT01, FFTW_REDFT01, FFTW_RODFT11, FFTW_REDFT11, FFTW_HC2R
    };

    fftw_r2r_kind kind_fwd_x[] = {transforms[BC.x()-1]};
    fftw_r2r_kind kind_fwd_y[] = {transforms[BC.y()-1]};
    fftw_r2r_kind kind_bwd_x[] = {transforms[BC.x()+4]};
    fftw_r2r_kind kind_bwd_y[] = {transforms[BC.y()+4]};

    // Input data containing random values

    int seed = 123 * rank;
    srand(seed);

    scalarBlock inputData(Ni[rank]);

    scalar average = 0;

    for ( int i = 0; i < Ni[rank].x(); i++ )
    {
        for ( int j = 0; j < Ni[rank].y(); j++ )
        {
            for ( int k = 0; k < Ni[rank].z(); k++ )
            {
                // inputData(i,j,k) =
                //   (si[rank].x() + i) * N.y() * N.z()
                // + (si[rank].y() + j) * N.z()
                // + (si[rank].z() + k);
                inputData(i,j,k) = rand()/RAND_MAX;
                average += inputData(i,j,k) / (cmptProduct(Ni[rank]));
            }
        }
    }

    for ( int i = 0; i < Ni[rank].x(); i++ )
    {
        for ( int j = 0; j < Ni[rank].y(); j++ )
        {
            for ( int k = 0; k < Ni[rank].z(); k++ )
            {
                inputData(i,j,k) -= average;
            }
        }
    }

    // Decompose into x-pencils

    labelVector X = decomp->X_;

    List<labelVector> Nx = decomp->Nx_;
    List<labelVector> sx = decomp->sx_;

    scalarBlock xData(Nx[rank]);

    // Transpose the data from the initial decomposition to x-pencils

    auto t1 = high_resolution_clock::now();

    decomp->transpose(inputData, xData, I, X);

    // FFT of xData in x-direction

    fft_rank = 1;
    howmany = Nx[rank].y() * Nx[rank].z();
    stride = Nx[rank].y() * Nx[rank].z();
    sep = 1;
    const int Na_x[] = {N.x()};

    plan_fwd_x = fftw_plan_many_r2r
    (
        fft_rank,
        Na_x,
        howmany,
        reinterpret_cast<double*>(xData.begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(xData.begin()),
        Na_x,
        stride,
        sep,
        kind_fwd_x,
        FFTW_ESTIMATE
    );

    fftw_execute(plan_fwd_x);

    // Decompose into y-pencils

    labelVector Y = decomp->Y_;

    List<labelVector> Ny = decomp->Ny_;
    List<labelVector> sy = decomp->sy_;

    scalarBlock yData(Ny[rank]);

    // Transpose the data

    decomp->transpose(xData, yData, X, Y);

    decomp->yTransFwd(yData);

    // FFT in y-direction

    fft_rank = 1;
    howmany = Ny[rank].x() * Ny[rank].z();
    stride = 1;
    sep = N.y();
    const int Na_y[] = {N.y()};

    plan_fwd_y = fftw_plan_many_r2r
    (
        fft_rank,
        Na_y,
        howmany,
        reinterpret_cast<double*>(yData.begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(yData.begin()),
        Na_y,
        stride,
        sep,
        kind_fwd_y,
        FFTW_ESTIMATE
    );

    fftw_execute(plan_fwd_y);

    decomp->yTransBwd(yData);

    // Decompose into z-pencils

    labelVector Z = decomp->Z_;

    List<labelVector> Nz = decomp->Nz_;
    List<labelVector> sz = decomp->sz_;

    scalarBlock zData(Nz[rank]);

    // Transpose the data

    decomp->transpose(yData, zData, Y, Z);

    // Eigenvalues lambda x

    scalar pi = 3.14159265359;

    List<scalar> lambda_x(Nz[rank].x());

    switch ( BC.x() )
    {
        case 1:
            for ( int i = 0; i < Nz[rank].x(); i++ )
            {
                lambda_x[i] = - 4.0
                * pow
                (
                    sin( (sz[rank].x()+i+1) * pi / (2.0 * N.x()) ),
                    2.0
                );
            }
            break;

        case 2:
            for ( int i = 0; i < Nz[rank].x(); i++ )
            {
                lambda_x[i] = - 4.0
                * pow
                (
                    sin( (sz[rank].x()+i) * pi / (2.0 * N.x()) ),
                    2.0
                );
            }
            break;

        case 3:
        case 4:
            for ( int i = 0; i < Nz[rank].x(); i++ )
            {
                lambda_x[i] = - 4.0
                * pow
                (
                    sin( (2.0*(sz[rank].x()+i+1)-1.0) * pi / (4.0 * N.x()) ),
                    2.0
                );
            }
            break;

        case 5:
            for ( int i = 0; i < Nz[rank].x(); i++ )
            {
                lambda_x[i] = -4.0
                * pow
                (
                    sin( (sz[rank].x()+i) * pi / N.x() ),
                    2.0
                );
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    // Eigenvalues lambda y

    List<scalar> lambda_y(Nz[rank].y());

    switch ( BC.y() )
    {
        case 1:
            for ( int j = 0; j < Nz[rank].y(); j++ )
            {
                lambda_y[j] = - 4.0
                * pow
                (
                    sin( (sz[rank].y()+j+1) * pi / (2.0 * N.y()) ),
                    2.0
                );
            }
            break;

        case 2:
            for ( int j = 0; j < Nz[rank].y(); j++ )
            {
                lambda_y[j] = - 4.0
                * pow
                (
                    sin( (sz[rank].y()+j) * pi / (2.0 * N.y()) ),
                    2.0
                );
            }
            break;

        case 3:
        case 4:
            for ( int j = 0; j < Nz[rank].y(); j++ )
            {
                lambda_y[j] = - 4.0
                * pow
                (
                    sin( (2.0*(sz[rank].y()+j+1)-1.0) * pi / (4.0 * N.y()) ),
                    2.0
                );
            }
            break;

        case 5:
            for ( int j = 0; j < Nz[rank].y(); j++ )
            {
                lambda_y[j] = -4.0
                * pow
                (
                    sin( (sz[rank].y()+j) * pi / N.y() ),
                    2.0
                );
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    // Solve for p_hat_hat

    scalarBlock pHat(Nz[rank]);

    List<scalar> D(N.z());

    switch ( BC.z() )
    {
        case 1:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    D[0]    = lambda_x[i] + lambda_y[j] - 3.0;
                    D[N.z()-1] = lambda_x[i] + lambda_y[j] - 3.0;
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i] + lambda_y[j] - 2.0;
                    }
                    tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(&zData(i,j,0)),
                        reinterpret_cast<double*>(&pHat(i,j,0))
                    );
                }
            }
            break;

        case 2:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    D[0]    = lambda_x[i] + lambda_y[j] - 1.0;
                    D[N.z()-1] = lambda_x[i] + lambda_y[j] - 1.0;
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i] + lambda_y[j] - 2.0;
                    }
                    // If coefficient matrix is singular
                    if ( lambda_x[i] == 0 && lambda_y[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(&zData(i,j,0)),
                        reinterpret_cast<double*>(&pHat(i,j,0))
                    );
                }
            }
            break;

        case 3:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    D[0]    = lambda_x[i] + lambda_y[j] - 3.0;
                    D[N.z()-1] = lambda_x[i] + lambda_y[j] - 1.0;
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i] + lambda_y[j] - 2.0;
                    }
                    tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(&zData(i,j,0)),
                        reinterpret_cast<double*>(&pHat(i,j,0))
                    );
                }
            }
            break;

        case 4:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    D[0]    = lambda_x[i] + lambda_y[j] - 1.0;
                    D[N.z()-1] = lambda_x[i] + lambda_y[j] - 3.0;
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i] + lambda_y[j] - 2.0;
                    }
                    tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(&zData(i,j,0)),
                        reinterpret_cast<double*>(&pHat(i,j,0))
                    );
                }
            }
            break;

        case 5:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    for (int k = 0; k < N.z(); k++ )
                    {
                        D[k] = lambda_x[i] + lambda_y[j] - 2.0;
                    }
                    // If coefficient matrix is singular
                    if ( lambda_x[i] == 0 && lambda_y[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    cyclicTridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(&zData(i,j,0)),
                        reinterpret_cast<double*>(&pHat(i,j,0))
                    );
                }
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    // Transpose to Y-pencils

    decomp->transpose(pHat, yData, Z, Y);

    // Inverse transform in y direction

    decomp->yTransFwd(yData);

    plan_bwd_y = fftw_plan_many_r2r
    (
        fft_rank,
        Na_y,
        howmany,
        reinterpret_cast<double*>(yData.begin()),
        Na_y,
        stride,
        sep,
        reinterpret_cast<double*>(yData.begin()),
        Na_y,
        stride,
        sep,
        kind_bwd_y,
        FFTW_ESTIMATE
    );

    fftw_execute(plan_bwd_y);

    decomp->yTransBwd(yData);

    // Transpose to x-pencils

    decomp->transpose(yData, xData, Y, X);

    // Inverse transform in x-direction

    fft_rank = 1;
    howmany = Nx[rank].y() * Nx[rank].z();
    stride = Nx[rank].y() * Nx[rank].z();
    sep = 1;

    plan_bwd_x = fftw_plan_many_r2r
    (
        fft_rank,
        Na_x,
        howmany,
        reinterpret_cast<double*>(xData.begin()),
        Na_x,
        stride,
        sep,
        reinterpret_cast<double*>(xData.begin()),
        Na_x,
        stride,
        sep,
        kind_bwd_x,
        FFTW_ESTIMATE
    );

    fftw_execute(plan_bwd_x);

    auto t2 = high_resolution_clock::now();

    xData /= normalization;

    // Gather p and f on main processor

    autoPtr<scalarBlock> pRecvBufferPtr;
    autoPtr<scalarBlock> fRecvBufferPtr;

    if ( ! rank )
    {
        pRecvBufferPtr.reset(new scalarBlock(N));
        fRecvBufferPtr.reset(new scalarBlock(N));
    }

    labelList recvCount(Pstream::nProcs());
    labelList recvDisplacement(Pstream::nProcs(), 0);

    // Gather p

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        recvCount[proc] = cmptProduct(Nx[proc])*sizeof(scalar);

        for (int p = 0; p < proc; p++)
            recvDisplacement[proc] += cmptProduct(Nx[p])*sizeof(scalar);
    }

    UPstream::gather
    (
        reinterpret_cast<char*>(xData.begin()),
        cmptProduct(Nx[rank])*sizeof(scalar),
        reinterpret_cast<char*>
        (
            rank == 0
          ? pRecvBufferPtr->begin()
          : nullptr
        ),
        recvCount,
        recvDisplacement,
        UPstream::worldComm
    );

    // Gather f

    for (label proc = 0; proc < Pstream::nProcs(); proc++)
    {
        recvCount[proc] = cmptProduct(Ni[proc])*sizeof(scalar);

        recvDisplacement[proc] = 0;

        for (int p = 0; p < proc; p++)
            recvDisplacement[proc] += cmptProduct(Ni[p])*sizeof(scalar);
    }

    UPstream::gather
    (
        reinterpret_cast<char*>(inputData.begin()),
        cmptProduct(Ni[rank])*sizeof(scalar),
        reinterpret_cast<char*>
        (
            rank == 0
          ? fRecvBufferPtr->begin()
          : nullptr
        ),
        recvCount,
        recvDisplacement,
        UPstream::worldComm
    );

    if ( ! rank )
    {
        scalarBlock p(N);
        scalarBlock f(N);

        decomp->unpack(pRecvBufferPtr(), sx, Nx, p, sx);
        decomp->unpack(fRecvBufferPtr(), si, Ni, f, si);

        if(check(p, f, BC))
        {
            Info << "-------------------------------------------" << nl;
            Info << "Pressure equation solution check successful" << nl;
            Info << "-------------------------------------------" << endl;

            Info<< "Completed in "
                << duration_cast<milliseconds>(t2 - t1).count()
                << " ms" << endl;
        }
    }
}
