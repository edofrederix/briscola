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

        const fvMesh& fvMsh = U.fvMsh();

        const colocatedVectorField& cc = fvMsh.metrics<colocated>().cellCenters();

        forAll(U, l)
        forAll(U[l], d)
        {
            const colocatedVectorDirection& ccd = cc[l][d];
            colocatedVectorDirection& Ud = U[l][d];

            forAllBlock(Ud, i, j, k)
            {
                Ud(i,j,k).x() =
                    ccd(i,j,k).y()
                + 0.4*Foam::sin(ccd(i,j,k).z()*3)
                + 0.1*Foam::sin(ccd(i,j,k).y()*3.1415927*2);

                Ud(i,j,k).z() =
                    Foam::sin(ccd(i,j,k).x()*3)
                * Foam::cos(ccd(i,j,k).y()*3.1415927);
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