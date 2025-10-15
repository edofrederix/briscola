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
    colocatedScalarField& p,
    colocatedScalarField& f,
    const fvMesh& fvMsh
)
{
    const FastPtrList<PartialList<scalar>>& cellSizes
        = fvMsh.msh().cast<rectilinearMesh>().globalCellSizes();

    labelVector Si =
        fvMsh.msh().decomp().globalStartPerProc()[Pstream::myProcNo()];

    PartialList<scalar> dx = cellSizes[0];
    PartialList<scalar> dy = cellSizes[1];
    PartialList<scalar> dz = cellSizes[2];

    forAllCells(p,i,j,k)
    {
        scalar dxw =
        (
            dx[Si.x() + i]
            + dx[Si.x() + i - 1]
        ) / 2.0;
        scalar dxe =
        (
            dx[Si.x() + i]
            + dx[Si.x() + i + 1]
        ) / 2.0;

        scalar dyw =
        (
            dy[Si.y() + j]
            + dy[Si.y() + j - 1]
        ) / 2.0;
        scalar dye =
        (
            dy[Si.y() + j]
            + dy[Si.y() + j + 1]
        ) / 2.0;

        scalar dzw =
        (
            dz[Si.z() + k]
            + dz[Si.z() + k - 1]
        ) / 2.0;
        scalar dze =
        (
            dz[Si.z() + k]
            + dz[Si.z() + k + 1]
        ) / 2.0;

        scalar residual =
              p(i-1,j,k) / (dxw*dx[Si.x() + i])
            - p(i,j,k)   / (dxw*dx[Si.x() + i])
            - p(i,j,k)   / (dxe*dx[Si.x() + i])
            + p(i+1,j,k) / (dxe*dx[Si.x() + i])

            + p(i,j-1,k) / (dyw*dy[Si.y() + j])
            - p(i,j,k)   / (dyw*dy[Si.y() + j])
            - p(i,j,k)   / (dye*dy[Si.y() + j])
            + p(i,j+1,k) / (dye*dy[Si.y() + j])

            + p(i,j,k-1) / (dzw*dz[Si.z() + k])
            - p(i,j,k)   / (dzw*dz[Si.z() + k])
            - p(i,j,k)   / (dze*dz[Si.z() + k])
            + p(i,j,k+1) / (dze*dz[Si.z() + k])

            + f(i,j,k);

        if(mag(residual) > 1e-10)
        {
            FatalError
                << "Test failed. Residual =  " << residual
                << " at index " << labelVector(i,j,k)
                << " on processor " << Pstream::myProcNo()
                << endl << abort(FatalError);
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
        true
    );

    colocatedScalarField p
    (
        "p",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true
    );

    labelVector N(fvMsh.msh().cast<rectilinearMesh>().N());

    Info << "Mesh size: " << N << endl;

    int seed = 123 * Pstream::myProcNo();
    srand(seed);

    scalar average = 0;

    forAllCells(f, i, j, k)
    {
        f(i,j,k) = 100.0 * static_cast<double>(rand()) / RAND_MAX - 0.5;
        average += f(i,j,k) / (cmptProduct(f.N()));
    }

    forAllCells(f, i, j, k)
    {
        f(i,j,k) -= average;
    }

    FFTPoissonSolver<stencil> solver(fvMsh);

    for (int r = 0; r < 1; r++)
    {
        solver.solve(p,f);
        Info << "Run number " << r+1 << " completed." << endl;
    }

    if(check(p, f, fvMsh))
    {
        Info << "-------------------------------------------" << nl;
        Info << "Pressure equation solution check successful" << nl;
        Info << "-------------------------------------------" << endl;
    }
}
