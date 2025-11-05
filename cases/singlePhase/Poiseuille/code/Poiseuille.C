#include "Poiseuille.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"

namespace Foam
{

using constant::mathematical::pi;

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(Poiseuille, 0);

addToRunTimeSelectionTable
(
    functionObject,
    Poiseuille,
    dictionary
);

bool Poiseuille::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        if (runTime_.foundObject<colocatedVectorField>("U"))
        {
            colocatedVectorField& U =
                runTime_.lookupObjectRef<colocatedVectorField>("U");

            const colocatedVectorField& cc =
                U.fvMsh().metrics<colocated>().cellCenters();

            forAllCells(U, i, j, k)
                U(i,j,k).x() = Foam::mag(cc(i,j,k).y()) < 1.0;

            U.correctBoundaryConditions();
        }
        else
        {
            staggeredScalarField& U =
                runTime_.lookupObjectRef<staggeredScalarField>("U");

            const fvMesh& fvMsh = U.fvMsh();

            const tensor base =
                fvMsh.msh().cast<rectilinearMesh>().base();

            const staggeredVectorField& cc =
                fvMsh.metrics<staggered>().cellCenters();

            forAllCells(U, d, i, j, k)
            {
                const vector u
                (
                    Foam::mag(cc(d,i,j,k).y()) < 1.0,
                    0,
                    0
                );

                U(d,i,j,k) = (base & u)[d];
            }

            U.correctBoundaryConditions();
        }
    }

    return true;
}

}

}

}

}