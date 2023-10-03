#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"
#include "constants.H"
#include "uniformMesh.H"
#include "exSchemes.H"
#include <fstream>

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

scalar initialCondition::TotalVolume_ = 0;
scalar initialCondition::BoundError_ = 0;

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

        colocatedVectorField& U =
            runTime_.lookupObjectRef<colocatedVectorField>("U");

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const fvMesh& fvMsh = U.fvMsh();

        const colocatedVectorField& cc =
            fvMsh.metrics<colocated>().cellCenters();

        const colocatedScalarField& cv =
            fvMsh.metrics<colocated>().cellVolumes();

        const meshField<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters();

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);

            U(i,j,k) = vector
                (
                  - 2.0*pi*(ccxy.y() - 0.5),
                    2.0*pi*(ccxy.x() - 0.5),
                    0
                );

            alpha(i,j,k) = 0;

            for (int l = 0; l < 8; l++)
            {
                vector v = vertex(i,j,k)[l];
                v.z() = 0;

                scalar tag = Foam::mag(v - vector(0.5,0.75,0)) < 0.15;

                if
                (
                    v.y() < 0.75-0.15+2*0.125
                 && v.x() > 0.5-0.025 && v.x() < 0.5+0.025
                )
                {
                    tag = 0;
                }

                alpha(i,j,k) += 0.125 * (tag);
            }

            TotalVolume_ += cv(i,j,k) * alpha(i,j,k);
        }

        reduce(TotalVolume_, sumOp<scalar>());

        alpha.correctBoundaryConditions();

        colocatedFaceScalarField& phi =
            runTime_.lookupObjectRef<colocatedFaceScalarField>("phi");

        U.correctBoundaryConditions();

        phi = ex::faceFlux(U);
    }

    return true;
}

bool initialCondition::execute()
{
    scalar LocalBoundError;

    colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const fvMesh& fvMsh = alpha.fvMsh();

    const colocatedScalarField& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    forAllCells(alpha, i, j, k)
    {
        LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);

        if (LocalBoundError > BoundError_)
            BoundError_ = LocalBoundError;
    }

    reduce(BoundError_, maxOp<scalar>());

    return true;
}

bool initialCondition::end()
{

    scalar L1Error = 0;
    scalar Volume = 0;
    scalar LocalBoundError;

    const colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const fvMesh& fvMsh = alpha.fvMsh();

    const colocatedScalarField& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    const meshField<vertexVector,colocated>& vertex =
        fvMsh.metrics<colocated>().vertexCenters();

    forAllCells(alpha, i, j, k)
    {
        scalar exactVol = 0;
        for (int l = 0; l < 8; l++)
        {
            vector v = vertex(i,j,k)[l];
            v.z() = 0;

            scalar tag = Foam::mag(v - vector(0.5,0.75,0)) < 0.15;

            if
            (
                v.y() < 0.75-0.15+2*0.125
             && v.x() > 0.5-0.025 && v.x() < 0.5+0.025
            )
            {
                tag = 0;
            }

            exactVol += 0.125 * (tag);
        }

        L1Error += cv(i,j,k) * Foam::mag(exactVol - alpha(i,j,k));
        Volume += cv(i,j,k) * alpha(i,j,k);
        LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);

        if (LocalBoundError > BoundError_)
            BoundError_ = LocalBoundError;
    }

    reduce(BoundError_, maxOp<scalar>());
    reduce(Volume, sumOp<scalar>());
    reduce(L1Error, sumOp<scalar>());

    Info<< "Total Volume: " << TotalVolume_ << nl
        << "Shape Error (absolute L1 norm): " << L1Error << nl
        << "Shape Error (relative L1 norm): " << L1Error / TotalVolume_ << nl
        << "Volume Error (absolute L1 norm): " << TotalVolume_ - Volume << nl
        << "Volume Error (relative L1 norm): "
        << (TotalVolume_ - Volume) / TotalVolume_ << nl
        << "Boundedness Error (infinite norm): " << BoundError_ << endl;

    return true;
}

}

}

}

}
