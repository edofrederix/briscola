#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"
#include "constants.H"
#include "exSchemes.H"

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
        const scalar pi = constant::mathematical::pi;

        colocatedVectorDirection& U =
            runTime_.lookupObjectRef<colocatedVectorField>("U")[0][0];

        colocatedScalarDirection& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

        const fvMesh& fvMsh = U.fvMsh();

        const colocatedVectorDirection& cc =
            fvMsh.metrics<colocated>().cellCenters()[0][0];

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);
            ccxy.z() = 0;

            U(i,j,k) = vector
                (
                  - 2.0*pi*(ccxy.y() - 0.5),
                    2.0*pi*(ccxy.x() - 0.5),
                    0
                );

            alpha(i,j,k) = Foam::mag(ccxy - vector(0.5,0.75,0)) < 0.15;

            if
            (
                ccxy.y() < 0.75-0.15+2*0.125
             && ccxy.x() > 0.5-0.025 && ccxy.x() < 0.5+0.025
            )
            {
                alpha(i,j,k) = 0;
            }
        }

        alpha.mshLevel().mshField().correctBoundaryConditions();

        colocatedFaceScalarField& phif =
            runTime_.lookupObjectRef<colocatedFaceScalarField>("phi");

        colocatedVectorField& Uf =
            runTime_.lookupObjectRef<colocatedVectorField>("U");

        Uf.correctBoundaryConditions();

        phif = ex::faceFlux(Uf);
    }

    return true;
}

}

}

}

}
