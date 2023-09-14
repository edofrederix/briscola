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

const scalar T = 3.0;
scalar TotalVolume = 0;
scalar BoundError = 0;

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

        const colocatedScalarDirection& cv =
            fvMsh.metrics<colocated>().cellVolumes()[0][0];

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);

            U(i,j,k) = vector
                (
                  2 * Foam::cos(pi * runTime_.time().value() / T) * Foam::sqr(Foam::sin(pi * ccxy.x())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * runTime_.time().value() / T) * Foam::sqr(Foam::sin(pi * ccxy.y())) * Foam::sin(2 * pi * ccxy.x()) * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * runTime_.time().value() / T) * Foam::sqr(Foam::sin(pi * ccxy.z())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.x())
                );

            alpha(i,j,k) = Foam::mag(ccxy - vector(0.35,0.35,0.35)) < 0.15;
            TotalVolume += cv(i,j,k) * alpha(i,j,k);
        }

        reduce(TotalVolume, sumOp<scalar>());

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

bool initialCondition::execute()
{
    const scalar pi = constant::mathematical::pi;
    scalar LocalBoundError;

    colocatedVectorDirection& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U")[0][0];

    colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const fvMesh& fvMsh = U.fvMsh();

    const colocatedVectorDirection& cc =
        fvMsh.metrics<colocated>().cellCenters()[0][0];

    const colocatedScalarDirection& cv =
        fvMsh.metrics<colocated>().cellVolumes()[0][0];

    forAllCells(U, i, j, k)
    {
        vector ccxy = cc(i,j,k);

        U(i,j,k) = vector
            (
                2 * Foam::cos(pi * (runTime_.time().value() + runTime_.deltaT().value())/ T) * Foam::sqr(Foam::sin(pi * ccxy.x())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.y())) * Foam::sin(2 * pi * ccxy.x()) * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.z())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.x())
            );

        LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);
        if (LocalBoundError > BoundError)
            BoundError = LocalBoundError;

    }

    reduce(BoundError, maxOp<scalar>());

    colocatedFaceScalarField& phif =
        runTime_.lookupObjectRef<colocatedFaceScalarField>("phi");

    colocatedVectorField& Uf =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    Uf.correctBoundaryConditions();

    phif = ex::faceFlux(Uf);

    return true;
}

bool initialCondition::end()
{
    if (Foam::mag(runTime_.time().value() - T) < 1e-12)
    {

        scalar L1Error = 0;
        scalar Volume = 0;
        scalar LocalBoundError;

        colocatedScalarDirection& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

        const fvMesh& fvMsh = alpha.fvMsh();

        const colocatedVectorDirection& cc =
            fvMsh.metrics<colocated>().cellCenters()[0][0];

        const colocatedScalarDirection& cv =
            fvMsh.metrics<colocated>().cellVolumes()[0][0];

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);

            L1Error += cv(i,j,k) * Foam::mag((Foam::mag(ccxy - vector(0.35,0.35,0.35)) < 0.15) - alpha(i,j,k));
            Volume += cv(i,j,k) * alpha(i,j,k);
            LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);

            if (LocalBoundError > BoundError)
                BoundError = LocalBoundError;

        }

        reduce(BoundError, maxOp<scalar>());
        reduce(Volume, sumOp<scalar>());
        reduce(L1Error, sumOp<scalar>());

        Info << "Total Volume: " << TotalVolume << endl <<
                "Shape Error (absolute L1 norm): " << L1Error << endl <<
                "Shape Error (relative L1 norm): " << L1Error / TotalVolume << endl <<
                "Volume Error (absolute L1 norm): " << TotalVolume - Volume << endl <<
                "Volume Error (relative L1 norm): " << (TotalVolume - Volume) / TotalVolume << endl <<
                "Boundness Error (infinite norm): " << BoundError << endl;
    }

    return true;
}

}

}

}

}
