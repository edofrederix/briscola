#include "channelFlow.H"
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

defineTypeNameAndDebug(channelFlow, 0);

addToRunTimeSelectionTable
(
    functionObject,
    channelFlow,
    dictionary
);

bool channelFlow::read(const dictionary& dict)
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
                    5300.0/360.0*2.0
                  * (
                        (1.0 - Foam::sqr(c.y()))
                      + 0.4*Foam::sin(c.z()*3)
                      + 0.1*Foam::sin(c.y()*pi)
                    );

                U(i,j,k).z() =
                    5300.0/360.0*2.0
                  * Foam::sin(c.x()*3)*Foam::sin(c.y()*pi);
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
                    5300.0/360.0*2.0
                  * (
                        (1.0 - Foam::sqr(c.y()))
                      + 0.4*Foam::sin(c.z()*3)
                      + 0.1*Foam::sin(c.y()*pi)
                    ),
                    0,
                    5300.0/360.0*2.0
                  * Foam::sin(c.x()*3)*Foam::sin(c.y()*pi)
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