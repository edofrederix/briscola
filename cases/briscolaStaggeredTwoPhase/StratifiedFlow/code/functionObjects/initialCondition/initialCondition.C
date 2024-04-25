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

        const meshField<vector,colocated>& cc =
            alpha.fvMsh().metrics<colocated>().cellCenters();

        const meshField<vector,staggered>& scc =
            U.fvMsh().metrics<staggered>().cellCenters();

        forAllCells(alpha, i, j, k)
        {
            // if (Foam::mag(cc(i,j,k)[1]) < 1e-3)
            // {
            //     alpha(i,j,k) = 0.5;
            // }
            // else
            if (cc(i,j,k)[1] < 0)
            {
                alpha(i,j,k) = 0;
            }
            else
            {
                alpha(i,j,k) = 1;
            }
        }

        forAllCells(U[0][0], i, j, k)
        {
            U(0,i,j,k) = 1.5 *
                        (
                            1.5 * (1 - Foam::sqr(scc(i,j,k)[1]/0.05))
                          + 0.25 * Foam::sin(4 * 3.14 * scc(i,j,k)[2])
                        );

            U(2,i,j,k) = 1.5 * 0.25 * Foam::sin(2 * 3.14 * scc(i,j,k)[0]);
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

