#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"
#include "constants.H"
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

        const meshDirection<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters()[0][0];

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

                scalar tag = Foam::mag(v - vector(0.5,0.75,0)) < 0.15;

                if
                (
                    v.y() < 0.75-0.15+0.125
                && v.x() > 0.5-0.025 && v.x() < 0.5+0.025
                )
                {
                    tag = 0;
                }
                alpha(i,j,k) += 0.125 * (tag);
            }

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
    scalar LocalBoundError;

    colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const fvMesh& fvMsh = alpha.fvMsh();

    const colocatedScalarDirection& cv =
        fvMsh.metrics<colocated>().cellVolumes()[0][0];

    forAllCells(alpha, i, j, k)
    {
        LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);
        if (LocalBoundError > BoundError)
            BoundError = LocalBoundError;

    }

    reduce(BoundError, maxOp<scalar>());

    return true;
}

bool initialCondition::end()
{

    scalar L1Error = 0;
    scalar Volume = 0;
    scalar LocalBoundError;

    colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const fvMesh& fvMsh = alpha.fvMsh();

    const colocatedScalarDirection& cv =
        fvMsh.metrics<colocated>().cellVolumes()[0][0];

    const meshDirection<vertexVector,colocated>& vertex =
        fvMsh.metrics<colocated>().vertexCenters()[0][0];

    forAllCells(alpha, i, j, k)
    {
        scalar exactVol = 0;

        for (int l = 0; l < 8; l++)
        {
            vector v = vertex(i,j,k)[l];

            scalar tag = Foam::mag(v - vector(0.5,0.75,0)) < 0.15;

            if
            (
                v.y() < 0.75-0.15+0.125
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

        if (LocalBoundError > BoundError)
            BoundError = LocalBoundError;

    }

    reduce(BoundError, maxOp<scalar>());
    reduce(Volume, sumOp<scalar>());
    reduce(L1Error, sumOp<scalar>());

    /*

    if (Pstream::master())
    {
        std::ofstream myfile;
        myfile.open ("errors.txt", std::ios::out|std::ios::app);
        myfile  << TotalVolume << " "
                << L1Error << " " << L1Error / TotalVolume << " "
                << TotalVolume - Volume << " " << (TotalVolume - Volume) / TotalVolume << " "
                << BoundError << "\n";
        myfile.close();
    }

    */

    Info << "Total Volume: " << TotalVolume << endl <<
            "Shape Error (absolute L1 norm): " << L1Error << endl <<
            "Shape Error (relative L1 norm): " << L1Error / TotalVolume << endl <<
            "Volume Error (absolute L1 norm): " << TotalVolume - Volume << endl <<
            "Volume Error (relative L1 norm): " << (TotalVolume - Volume) / TotalVolume << endl <<
            "Boundness Error (infinite norm): " << BoundError << endl;

    return true;
}

}

}

}

}
