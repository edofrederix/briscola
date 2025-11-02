#include "stratifiedFlow.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include <fstream>

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

defineTypeNameAndDebug(stratifiedFlow, 0);

addToRunTimeSelectionTable
(
    functionObject,
    stratifiedFlow,
    dictionary
);

bool stratifiedFlow::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        // Volume fraction

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const meshField<vector,colocated>& cc =
            alpha.fvMsh().metrics<colocated>().cellCenters();

        alpha = Zero;

        forAllCells(alpha, i, j, k)
            alpha(i,j,k) = cc(i,j,k).y() > 0;

        alpha.correctBoundaryConditions();

        // Velocity

        if (runTime_.foundObject<colocatedVectorField>("U"))
        {
            colocatedVectorField& U =
                runTime_.lookupObjectRef<colocatedVectorField>("U");

            U = Zero;

            forAllCells(U, i, j, k)
            {
                const vector c = cc(i,j,k);

                U(i,j,k).x() =
                    1.5
                  * (
                        1.5*(1.0 - Foam::sqr(c.y()/0.05))
                      + 0.25*Foam::sin(4.0*pi*c.z()/0.2)
                    );

                U(i,j,k).z() =
                    1.5*0.25*Foam::sin(2.0*pi*c.x()/0.4);
            }

            U.correctBoundaryConditions();
        }
        else
        {
            staggeredScalarField& U =
                runTime_.lookupObjectRef<staggeredScalarField>("U");

            U = Zero;

            const meshField<vector,staggered>& scc =
                alpha.fvMsh().metrics<staggered>().cellCenters();

            const tensor base =
                U.fvMsh().msh().cast<rectilinearMesh>().base();

            forAllCells(U, d, i, j, k)
            {
                const vector c = scc(d,i,j,k);

                const vector u
                (
                    1.5
                  * (
                        1.5*(1.0 - Foam::sqr(c.y()/0.05))
                      + 0.25*Foam::sin(4.0*pi*c.z()/0.2)
                    ),
                    0,
                    1.5*0.25*Foam::sin(2.0*pi*c.x()/0.4)
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

