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

        const fvMesh& fvMsh = U.fvMsh();

        const staggeredVectorField& cc =
            fvMsh.metrics<staggered>().cellCenters();

        forAll(U, l)
        {
            const staggeredVectorDirection& cc0 = cc[l][0];
            const staggeredVectorDirection& cc2 = cc[l][2];

            staggeredScalarDirection& U0 = U[l][0];
            staggeredScalarDirection& U2 = U[l][2];

            forAllCells(U0, i, j, k)
            {
                U0(i,j,k) =
                    5300.0/360.0*2.0
                  * (
                        (1.0-Foam::sqr(cc0(i,j,k).y()))
                      + 0.4*Foam::sin(cc0(i,j,k).z()*3)
                      + 0.1*Foam::sin(cc0(i,j,k).y()*3.1415927)
                    );
            }

            forAllCells(U2, i, j, k)
            {
                U2(i,j,k) =
                    5300.0/360.0/2.0
                  * Foam::sin(cc2(i,j,k).x()*3)
                  * Foam::sin(cc2(i,j,k).y()*3.1415927);
            }
        }

        U.correctBoundaryConditions();
    }

    return true;
}

}

}

}

}