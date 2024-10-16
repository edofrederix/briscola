#include "arguments.H"
#include "Time.H"
#include "FFTPoissonSolver.H"
#include "fv.H"
#include "rectilinearMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

bool check
(
    colocatedScalarField& x,
    colocatedScalarField& b,
    const fvMesh& fvMsh,
    const scalar deltaT
)
{
    labelVector Nf(x.B().shape());
    labelVector N(x.fvMsh().msh().cast<rectilinearMesh>().N());

    const PtrList<PartialList<scalar>>& cellSizes
        = fvMsh.msh().cast<rectilinearMesh>().globalCellSizes();

    labelVector Si =
        fvMsh.msh().decomp().globalStartPerProc()[Pstream::myProcNo()];

    scalarList dx2 = sqr(cellSizes[0]);
    scalarList dy2 = sqr(cellSizes[1]);
    scalarList dz2 = sqr(cellSizes[2]);

    for (int i = 1; i < Nf.x() - 1; i++)
    {
        for (int j = 1; j < Nf.y() - 1; j++)
        {
            for (int k = 1; k < Nf.z() - 1; k++)
            {
                scalar residual =
                (
                    x.B()(i-1,j,k)
                  - 2.0 * x.B()(i,j,k)
                  + x.B()(i+1,j,k)
                ) / dx2[Si.x() + i-1]
              + (
                    x.B()(i,j-1,k)
                  - 2.0 * x.B()(i,j,k)
                  + x.B()(i,j+1,k)
                ) / dy2[Si.y() + j-1]
              + (
                    x.B()(i,j,k-1)
                  - 2.0 * x.B()(i,j,k)
                  + x.B()(i,j,k+1)
                ) / dz2[Si.z() + k-1]
              - x.B()(i,j,k) / deltaT
              + b.B()(i,j,k)
              + x.oldTime().B()(i,j,k) / deltaT;

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

    colocatedScalarField b
    (
        "b",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    colocatedScalarField x
    (
        "x",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    labelVector N(fvMsh.msh().cast<rectilinearMesh>().N());

    Info << "Mesh size: " << N << endl;

    int seed = 123 * Pstream::myProcNo();
    srand(seed);

    forAllCells(b, i, j, k)
    {
        b(i,j,k) = static_cast<double>(rand()) / RAND_MAX - 0.5;
        x(i,j,k) = static_cast<double>(rand()) / RAND_MAX - 0.5;
    }

    x.setOldTime();

    FFTPoissonSolver<stencil> solver(fvMsh);

    solver.solve(x, b, true);

    const scalar deltaT = runTime.deltaTValue();

    if(check(x, b, fvMsh, deltaT))
    {
        Info << "-------------------------------------------" << nl;
        Info << "Pressure equation solution check successful" << nl;
        Info << "-------------------------------------------" << endl;
    }
}
