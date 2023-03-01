#include "FFTPoissonSolver.H"
#include "imSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void FFTPoissonSolver::checkMesh() const
{
    if (this->fvMsh_.rectilinear() != unitXYZ)
    {
        FatalErrorInFunction
            << "Mesh must be rectilinear in three directions "
            << "for the " << this->type() << " solver." << endl
            << abort(FatalError);
    }

    if (cmptSum(this->fvMsh_.uniform()) < 2)
    {
        FatalErrorInFunction
            << "At least two mesh directions must be uniform "
            << "for the " << this->type() << " solver." << endl
            << abort(FatalError);
    }
}

void FFTPoissonSolver::prepare()
{
    cellSizes_.clear();
    cellSizes_.setSize(3);

    for (int dir = 0; dir < 3; dir++)
    {
        cellSizes_.set
        (
            dir,
            new scalarList(this->fvMsh_.rectilinearCellSizes(dir))
        );
    }
}

FFTPoissonSolver::FFTPoissonSolver
(
    const word PoissonSolverName,
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dict,fvMsh)
{
    checkMesh();
    prepare();

    decomp_ = new decomposer(fvMsh);
    tds_ = new tridiagonalSolver(fvMsh.N().z());
}

FFTPoissonSolver::FFTPoissonSolver
(
    const fvMesh& fvMsh
)
:
    PoissonSolver<stencil,scalar,colocated>(dictionary(),fvMsh)
{
    checkMesh();
    prepare();

    decomp_ = new decomposer(fvMsh);
    tds_ = new tridiagonalSolver(fvMsh.N().z());
}

void FFTPoissonSolver::solve
(
    colocatedScalarField& x,
    const colocatedScalarField* bPtr,
    const colocatedFaceScalarField* lambdaPtr,
    const bool ddt
)
{
    using Foam::pow;
    using Foam::sin;

    // Mesh dimensions

    labelVector N(this->fvMsh_.N());

    // Global boundary conditions
    // TO-DO: Change the whole confusing BC vector concept

    labelVector BC = vector(-1,-1,-1);

    switch (globalBoundaryConditionBaseType(x, faceOffsets[0]))
    {
        case 3:
            BC.x() = 5; // C-C
            break;
        case 4:
            if (globalBoundaryConditionBaseType(x, faceOffsets[1]) == 4)
            {
                BC.x() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[1]) == 5)
            {
                BC.x() = 3; // D-N
            }
            break;
        case 5:
            if (globalBoundaryConditionBaseType(x, faceOffsets[1]) == 4)
            {
                BC.x() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[1]) == 5)
            {
                BC.x() = 2; // N-N
            }
            break;
        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(x, faceOffsets[2]))
    {
        case 3:
            BC.y() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(x, faceOffsets[3]) == 4)
            {
                BC.y() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[3]) == 5)
            {
                BC.y() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(x, faceOffsets[3]) == 4)
            {
                BC.y() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[3]) == 5)
            {
                BC.y() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }

    switch (globalBoundaryConditionBaseType(x, faceOffsets[4]))
    {
        case 3:
            BC.z() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(x, faceOffsets[5]) == 4)
            {
                BC.z() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[5]) == 5)
            {
                BC.z() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(x, faceOffsets[5]) == 4)
            {
                BC.z() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(x, faceOffsets[5]) == 5)
            {
                BC.z() = 2; // N-N
            }
            break;

        default:
            FatalError
                << "Incorrect pressure boundary condition." << endl
                << abort(FatalError);
            break;
    }

    // FFT variables

    int fft_rank, howmany, stride, sep;

    fftw_plan plan_fwd_x;
    fftw_plan plan_fwd_y;

    fftw_plan plan_bwd_x;
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

    // Initial decomposition

    labelVector I(decomp_->I()); // TO-DO: This is only for one brick

    List<labelVector> Ni = decomp_->Ni();
    List<labelVector> si = decomp_->si();

    // Decompose into x-pencils

    labelVector X = decomp_->X();

    List<labelVector> Nx = decomp_->Nx();
    List<labelVector> sx = decomp_->sx();

    int rank = Pstream::myProcNo();
    scalarBlock xData(Nx[rank]);

    // Copy the data from bPtr (the RHS) to a scalarBlock

    scalarBlock initData(Ni[rank]);

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                initData(i,j,k) = bPtr[0][0][0](i,j,k);
            }
        }
    }

    // Transpose the RHS to x-pencils


    decomp_->transpose(initData, xData, I, X);

    // RHS = - bPtr

    xData *= -1.0;

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

    labelVector Y = decomp_->Y();

    List<labelVector> Ny = decomp_->Ny();
    List<labelVector> sy = decomp_->sy();

    scalarBlock yData(Ny[rank]);

    // Transpose the data

    decomp_->transpose(xData, yData, X, Y);

    decomp_->yTransFwd(yData);

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

    decomp_->yTransBwd(yData);

    // Decompose into z-pencils

    labelVector Z = decomp_->Z();

    List<labelVector> Nz = decomp_->Nz();
    List<labelVector> sz = decomp_->sz();

    scalarBlock zData(Nz[rank]);

    // Transpose the data

    decomp_->transpose(yData, zData, Y, Z);

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

    scalarList dx2 = sqr(cellSizes_[0]);
    scalarList dy2 = sqr(cellSizes_[1]);
    scalarList dz2 = sqr(cellSizes_[2]);

    List<scalar> DU = 1.0/dz2;
    List<scalar> D(N.z());
    List<scalar> DL = 1.0/dz2;

    switch ( BC.z() )
    {
        case 1:
            for (int i = 0; i < Nz[rank].x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz[rank].y(); j++) // solve column by column
                {
                    D[0]       = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 3.0/dz2[0]; // TO-DO nonuniform mesh?
                    D[N.z()-1] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 3.0/dz2[0];
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 2.0/dz2[0];
                    }
                    tds_->tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
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
                    D[0]       = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 1.0/dz2[0];
                    D[N.z()-1] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 1.0/dz2[0];
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 2.0/dz2[0];
                    }
                    // If coefficient matrix is singular
                    if ( lambda_x[i] == 0 && lambda_y[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    tds_->tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
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
                    D[0]       = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 3.0/dz2[0];
                    D[N.z()-1] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 1.0/dz2[0];
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 2.0/dz2[0];
                    }
                    tds_->tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
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
                    D[0]       = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 1.0/dz2[0];
                    D[N.z()-1] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 3.0/dz2[0];
                    for (int k = 1; k < N.z()-1; k++ )
                    {
                        D[k] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 2.0/dz2[0];
                    }
                    tds_->tridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
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
                        D[k] = lambda_x[i]/dx2[0] + lambda_y[j]/dy2[0] - 2.0/dz2[0];
                    }
                    // If coefficient matrix is singular
                    if ( lambda_x[i] == 0 && lambda_y[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    tds_->cyclicTridiagonalSolve
                    (
                        N.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
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

    decomp_->transpose(pHat, yData, Z, Y);

    // Inverse transform in y direction

    decomp_->yTransFwd(yData);

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

    decomp_->yTransBwd(yData);

    // Transpose to x-pencils

    decomp_->transpose(yData, xData, Y, X);

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

    xData /= normalization;

    // Transpose back to the p-field

    scalarBlock pData(Ni[rank]);

    decomp_->transpose(xData, pData, X, I);

    // Copy the p scalarBlock to the x meshField

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                x[0][0](i,j,k) = pData(i,j,k);
            }
        }
    }

    x.correctBoundaryConditions();

    Info << "FFT: Solving for colocated p, residual = 0, nIter = 1" << endl;
}

}

}

}
