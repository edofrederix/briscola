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
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const colocatedVectorField& cc =
            alpha.fvMsh().metrics<colocated>().cellCenters();

        forAllCells(alpha, i, j, k)
            if
            (
                Foam::sqr(cc(i,j,k).x()) + Foam::sqr(cc(i,j,k).y()-0.5)
              < Foam::sqr(0.25)
            )
                alpha(i,j,k) = 1.0;

        alpha.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}
