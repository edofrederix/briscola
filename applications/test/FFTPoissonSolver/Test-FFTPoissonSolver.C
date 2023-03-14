#include "arguments.H"
#include "Time.H"
#include "FFTPoissonSolver.H"
#include "fv.H"
#include "rectilinearMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

// Check the pressure field including ghost cells

bool check(scalarBlock& p, scalarBlock& f, labelVector BC)
{
    labelVector N(p.shape());

    // uniform meshes only
    scalar dx2 = sqr(1.0/N.x());
    scalar dy2 = sqr(1.0/N.y());
    scalar dz2 = sqr(1.0/N.z());


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
                (
                  ((i == 0) ?          px0 : p(i-1,j,k))
                  - 2.0 * p(i,j,k)
                  + ((i == N.x() - 1) ?  pxN : p(i+1,j,k))
                ) / dx2
                + (
                  ((j == 0) ?          py0 : p(i,j-1,k))
                  - 2.0 * p(i,j,k)
                  + ((j == N.y() - 1) ?  pyN : p(i,j+1,k))
                ) / dy2
                + (
                  ((k == 0) ?          pz0 : p(i,j,k-1))
                  - 2.0 * p(i,j,k)
                  + ((k == N.z() - 1) ?  pzN : p(i,j,k+1))
                ) / dz2
                - f(i,j,k);

                if(mag(residual) > 1e-10)
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

bool check(colocatedScalarField& p, colocatedScalarField& f)
{
    labelVector Nf(p[0][0].B().shape());
    labelVector N(p.fvMsh().msh().cast<rectilinearMesh>().N());

    // uniform 1m x 1m x 1m meshes only
    scalar dx2 = sqr(1.0/N.x());
    scalar dy2 = sqr(1.0/N.y());
    scalar dz2 = sqr(1.0/N.z());


    // scalar f;

    for (int i = 1; i < Nf.x() - 1; i++)
    {
        for (int j = 1; j < Nf.y() - 1; j++)
        {
            for (int k = 1; k < Nf.z() - 1; k++)
            {
                scalar residual =
                (
                    p[0][0].B()(i-1,j,k) - 2.0 * p[0][0].B()(i,j,k) + p[0][0].B()(i+1,j,k)
                ) / dx2
                + (
                    p[0][0].B()(i,j-1,k) - 2.0 * p[0][0].B()(i,j,k) + p[0][0].B()(i,j+1,k)
                ) / dy2
                + (
                    p[0][0].B()(i,j,k-1) - 2.0 * p[0][0].B()(i,j,k) + p[0][0].B()(i,j,k+1)
                ) / dz2
                + f[0][0].B()(i,j,k);

                if(mag(residual) > 1e-10)
                {
                    FatalError
                        << "Test failed. Residual =  " << residual
                        << " at index " << labelVector(i,j,k)
                        << " on processor " << Pstream::myProcNo()
                        << endl << abort(FatalError);
                }
            }
        }
    }

    return true;
}


int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    using std::rand;
    using std::srand;

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    colocatedScalarField f
    (
        "f",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    );

    colocatedScalarField p
    (
        "p",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    );

    for (int i = 0; i < 6; i++)
    {
        const labelVector bo = faceOffsets[i];

        Info << "Global BC " << i << bo << ": " <<globalBoundaryConditionBaseType(p, bo) << endl;
    }

    labelVector BC(-1,-1,-1);

    switch (globalBoundaryConditionBaseType(p, faceOffsets[0]))
    {
        case 3:
            BC.x() = 5; // C-C
            break;
        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 4)
            {
                BC.x() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 5)
            {
                BC.x() = 3; // D-N
            }
            break;
        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 4)
            {
                BC.x() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[1]) == 5)
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

    switch (globalBoundaryConditionBaseType(p, faceOffsets[2]))
    {
        case 3:
            BC.y() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 4)
            {
                BC.y() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 5)
            {
                BC.y() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 4)
            {
                BC.y() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[3]) == 5)
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

    switch (globalBoundaryConditionBaseType(p, faceOffsets[4]))
    {
        case 3:
            BC.z() = 5; // C-C
            break;

        case 4:
            if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 4)
            {
                BC.z() = 1; // D-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 5)
            {
                BC.z() = 3; // D-N
            }
            break;

        case 5:
            if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 4)
            {
                BC.z() = 4; // N-D
            }
            else if (globalBoundaryConditionBaseType(p, faceOffsets[5]) == 5)
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

    Info << "BC vector: " << BC << endl;

    autoPtr<decomposer> decomp;
    decomp = new decomposer(fvMsh);

    labelVector N(fvMsh.msh().cast<rectilinearMesh>().N());

    Info << "Mesh size: " << N << endl;

    labelVector I(decomp->I());

    List<labelVector> Ni = decomp->Ni();
    List<labelVector> si = decomp->si();

    Info << "Initial decomposition: " << I << endl;

    int seed = 123 * Pstream::myProcNo();
    srand(seed);

    int rank = Pstream::myProcNo();

    scalar average = 0;

    scalarBlock fCopy(Ni[rank]);

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                f[0][0](i,j,k) = 100.0 * static_cast<double>(rand()) / RAND_MAX - 0.5;
                average += f[0][0](i,j,k) / (cmptProduct(Ni[rank]));
            }
        }
    }

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                f[0][0](i,j,k) -= average;
                fCopy(i,j,k) = - f[0][0](i,j,k); // minus sign because 0 = laplacian(p) - f
            }
        }
    }

    FFTPoissonSolver solver(fvMsh);

    Info << "FFT Solver constructed. " << endl;

    for (int r = 0; r < 1; r++)
    {
        solver.solve(p,f);
        Info << "Run number " << r+1 << " completed." << endl;
    }

    scalarBlock pCopy(Ni[rank]);

    for (int i = 0; i < Ni[rank].x(); i++)
    {
        for (int j = 0; j < Ni[rank].y(); j++)
        {
            for (int k = 0; k < Ni[rank].z(); k++)
            {
                pCopy(i,j,k) = p[0][0](i,j,k);
            }
        }
    }

    // for (int i = 0; i < p[0][0].B().shape().x(); i++)
    // {
    //     for (int j = 0; j < p[0][0].B().shape().y(); j++)
    //     {
    //         for (int k = 0; k < p[0][0].B().shape().z(); k++)
    //         {
    //             Info << p[0][0].B()(i,j,k) << "     ";
    //         }
    //         Info << nl;
    //     }
    //     Info << nl;
    // }
    // Info << endl;

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
        recvCount[proc] = cmptProduct(Ni[proc])*sizeof(scalar);

        for (int pr = 0; pr < proc; pr++)
            recvDisplacement[proc] += cmptProduct(Ni[pr])*sizeof(scalar);
    }

    UPstream::gather
    (
        reinterpret_cast<char*>(pCopy.begin()),
        cmptProduct(Ni[rank])*sizeof(scalar),
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

        for (int pr = 0; pr < proc; pr++)
            recvDisplacement[proc] += cmptProduct(Ni[pr])*sizeof(scalar);
    }

    UPstream::gather
    (
        reinterpret_cast<char*>(fCopy.begin()),
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
        scalarBlock pFull(N);
        scalarBlock fFull(N);

        decomp->unpack(pRecvBufferPtr(), si, Ni, pFull, si);
        decomp->unpack(fRecvBufferPtr(), si, Ni, fFull, si);

        if(check(pFull, fFull, BC))
        {
            Info << "-------------------------------------------" << nl;
            Info << "Pressure equation solution check successful" << nl;
            Info << "-------------------------------------------" << endl;
        }
    }

    if(check(p, f))
        {
            Info << "-------------------------------------------" << nl;
            Info << "Pressure equation solution check successful" << nl;
            Info << "-------------------------------------------" << endl;
        }
}
