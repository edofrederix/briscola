#include "pipeFlow.H"
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

defineTypeNameAndDebug(pipeFlow, 0);

addToRunTimeSelectionTable
(
    functionObject,
    pipeFlow,
    dictionary
);

bool pipeFlow::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedVectorField& U =
            runTime_.lookupObjectRef<colocatedVectorField>("U");

        const colocatedVectorField& cc =
            U.fvMsh().metrics<colocated>().cellCenters();

        forAllCells(U, i, j, k)
        {
            const scalar x = cc(i,j,k).x();
            const scalar y = cc(i,j,k).y();
            const scalar z = cc(i,j,k).z();

            const scalar r = Foam::sqrt(Foam::sqr(x)+Foam::sqr(y));
            const scalar theta = Foam::atan2(x,y);

            U(i,j,k).x() =
                5300.0/360.0*Foam::sin(r*pi*2.0)
              * Foam::sin(z*2.0)*Foam::cos(theta)/4.0;

            U(i,j,k).y() =
                5300.0/360.0*Foam::sin(r*pi*2.0)
              * Foam::sin(z*2.0)*Foam::sin(theta)/4.0;

            U(i,j,k).z() =
                5300.0/360.0
              * (3.0*(1.0-Foam::sqr(r)) + Foam::sin(r*8.0*pi)/4.0);
        }

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
