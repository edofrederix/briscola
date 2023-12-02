#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(initialCondition, 0);

addToRunTimeSelectionTable
(
    functionObject,
    initialCondition,
    dictionary
);

initialCondition::initialCondition
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{
    read(dict);
}

initialCondition::~initialCondition()
{}

bool initialCondition::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedVectorField& U =
            runTime_.lookupObjectRef<colocatedVectorField>("U");

        const colocatedVectorField& cc =
            U.fvMsh().metrics<colocated>().cellCenters();

        const scalar pi = 3.1415927;

        forAllLevels(U, l, d, i, j, k)
        {
            const scalar x = cc(l,d,i,j,k).x();
            const scalar y = cc(l,d,i,j,k).y();
            const scalar z = cc(l,d,i,j,k).z();

            const scalar r = Foam::sqrt(Foam::sqr(x)+Foam::sqr(y));
            const scalar theta = Foam::atan2(x,y);

            U(l,d,i,j,k).x() =
                5300.0/360.0*Foam::sin(r*pi*2.0)
              * Foam::sin(z*2.0)*Foam::cos(theta)/4.0;

            U(l,d,i,j,k).y() =
                5300.0/360.0*Foam::sin(r*pi*2.0)
              * Foam::sin(z*2.0)*Foam::sin(theta)/4.0;

            U(l,d,i,j,k).z() =
                5300.0/360.0
              * (2.0*(1.0-Foam::sqr(r)) + Foam::sin(r*4.0*pi)/4.0);
        }

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
