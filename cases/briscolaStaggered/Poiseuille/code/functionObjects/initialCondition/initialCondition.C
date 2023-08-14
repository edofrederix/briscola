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

        U = Zero;

        forAll(U, l)
            U[l][0] = 1.0;

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
