#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include <fstream>

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

        staggeredScalarField& U =
            runTime_.lookupObjectRef<staggeredScalarField>("U");

        alpha = Zero;
        U = Zero;

        const meshField<vector,colocated>& ccc =
            alpha.fvMsh().metrics<colocated>().cellCenters();

        const meshField<vector,staggered>& scc =
            alpha.fvMsh().metrics<staggered>().cellCenters();

        forAllCells(alpha, l, d, i, j, k)
            alpha(l,d,i,j,k) = ccc(l,d,i,j,k).y() < 0.0 ? 0.0 : 1.0;

        forAllCells(U, l, d, i, j, k)
        {
            if (d == 0)
                U(l,d,i,j,k) =
                    1.5
                  * (
                        1.5*(1.0 - Foam::sqr(scc(l,d,i,j,k).y()/0.05))
                      + 0.25*Foam::sin(4.0*3.14*scc(i,j,k).z())
                    );

            if (d == 2)
                U(l,d,i,j,k) = 1.5*0.25*Foam::sin(2.0*3.14*scc(l,d,i,j,k).x());
        }

        alpha.correctBoundaryConditions();
        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}

