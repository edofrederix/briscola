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

        forAllCells(U, l, d, i, j, k)
        {
            U(l,d,i,j,k).x() =
                cc(l,d,i,j,k).y()
              + 0.4*Foam::sin(cc(l,d,i,j,k).z()*3)
              + 0.1*Foam::sin(cc(l,d,i,j,k).y()*3.1415927*2);

            U(l,d,i,j,k).z() =
                Foam::sin(cc(l,d,i,j,k).x()*3)
              * Foam::cos(cc(l,d,i,j,k).y()*3.1415927);
        }

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}