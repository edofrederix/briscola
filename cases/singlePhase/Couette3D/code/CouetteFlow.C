#include "CouetteFlow.H"
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

defineTypeNameAndDebug(CouetteFlow, 0);

addToRunTimeSelectionTable
(
    functionObject,
    CouetteFlow,
    dictionary
);

bool CouetteFlow::read(const dictionary& dict)
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
            {
                const vector c = cc(i,j,k);

                U(i,j,k).x() =
                    cc(i,j,k).y()
                  + 0.4*Foam::sin(c.z()*3)
                  + 0.1*Foam::sin(c.y()*pi*2);

                U(i,j,k).z() =
                    Foam::sin(c.x()*3)*Foam::cos(c.y()*pi);
            }

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
                const vector c = cc(d,i,j,k);

                const vector u
                (
                    cc(i,j,k).y()
                  + 0.4*Foam::sin(c.z()*3)
                  + 0.1*Foam::sin(c.y()*pi*2),
                    0,
                    Foam::sin(c.x()*3)*Foam::cos(c.y()*pi)
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