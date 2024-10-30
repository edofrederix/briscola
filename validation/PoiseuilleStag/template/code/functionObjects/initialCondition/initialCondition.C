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
        staggeredScalarField& U =
            runTime_.lookupObjectRef<staggeredScalarField>("U");

        const colocatedVectorField& cc =
            U.fvMsh().metrics<colocated>().cellCenters();

        forAllCells(U, l, d, i, j, k)
            if (d == 0)
                U(l,d,i,j,k) = Foam::mag(cc(l,d,i,j,k).y()) < 1.0;

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
