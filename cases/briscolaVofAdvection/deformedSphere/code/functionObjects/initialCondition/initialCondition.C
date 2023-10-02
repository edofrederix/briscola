#include "initialCondition.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"
#include "uniformMesh.H"
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

const scalar T = 3.0;
scalar TotalVolume = 0;
scalar BoundError = 0;

int Nx = 10;
int Ny = 10;
int Nz = 10;

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

        const mesh& mh = fvMsh.msh();

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);

            U(i,j,k) = vector
                (
                  2 * Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.x())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.y())) * Foam::sin(2 * pi * ccxy.x()) * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.z())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.x())
                );
        }

        if (mh.castable<uniformMesh>())
        {
            const uniformMesh& uniMh = mh.cast<uniformMesh>();
            vector cellSizes = uniMh.cellSize();

            vertexVector refCell;

            refCell.lba() = vector(0, 0, 0);
            refCell.rba() = vector(cellSizes[0]/double(Nx), 0, 0);
            refCell.lta() = vector(0, cellSizes[1]/double(Ny), 0);
            refCell.rta() = vector(cellSizes[0]/double(Nx), cellSizes[1]/double(Ny), 0);
            refCell.lbf() = vector(0, 0, cellSizes[2]/double(Nz));
            refCell.rbf() = vector(cellSizes[0]/double(Nx), 0, cellSizes[2]/double(Nz));
            refCell.ltf() = vector(0, cellSizes[1]/double(Ny), cellSizes[2]/double(Nz));
            refCell.rtf() = vector(cellSizes[0]/double(Nx), cellSizes[1]/double(Ny), cellSizes[2]/double(Nz));

            forAllCells(alpha, i, j, k)
            {

                alpha(i,j,k) = 0;

                for (int aux1 = 0; aux1 < Nx; aux1++)
                {
                    for (int aux2 = 0; aux2 < Ny; aux2++)
                    {
                        for (int aux3 = 0; aux3 < Nz; aux3++)
                        {
                            vertexVector cell = vertex(i,j,k).lba()
                                + (double(aux1)) * vector(cellSizes[0]/double(Nx), 0, 0)
                                + (double(aux2)) * vector(0, cellSizes[1]/double(Ny), 0)
                                + (double(aux3)) * vector(0, 0, cellSizes[2]/double(Nz))
                                + refCell;

                            for (int l = 0; l < 8; l++)
                            {
                                vector v = cell[l];
                                alpha(i,j,k) += (1/(double(Nx * Ny * Nz) * cellSizes[0] * cellSizes[1] * cellSizes[2])) * 0.125
                                        * Foam::min(cellSizes[0] * cellSizes[1] * cellSizes[2], Foam::max(0,
                                        10 * (Foam::sqr(0.15) - Foam::sqr(Foam::mag(v - vector(0.35,0.35,0.35))))));
                            }
                        }
                    }
                }

                alpha(i,j,k) = Foam::min(1, alpha(i,j,k));
                alpha(i,j,k) = Foam::max(0, alpha(i,j,k));

                TotalVolume += cv(i,j,k) * alpha(i,j,k);
            }
        }
        else
        {
            forAllCells(alpha, i, j, k)
            {
                alpha(i,j,k) = 0;
                for (int l = 0; l < 8; l++)
                {
                    vector v = vertex(i,j,k)[l];
                    alpha(i,j,k) += 0.125 * (Foam::mag(v - vector(0.35,0.35,0.35)) < 0.15);
                }
                TotalVolume += cv(i,j,k) * alpha(i,j,k);
            }
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
                2 * Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value())/ T) * Foam::sqr(Foam::sin(pi * ccxy.x())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.y())) * Foam::sin(2 * pi * ccxy.x()) * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T) * Foam::sqr(Foam::sin(pi * ccxy.z())) * Foam::sin(2 * pi * ccxy.y()) * Foam::sin(2 * pi * ccxy.x())
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

        const colocatedScalarDirection& cv =
            fvMsh.metrics<colocated>().cellVolumes()[0][0];

        const meshDirection<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters()[0][0];

        const mesh& mh = fvMsh.msh();

        if (mh.castable<uniformMesh>())
        {
            const uniformMesh& uniMh = mh.cast<uniformMesh>();
            vector cellSizes = uniMh.cellSize();

            vertexVector refCell;

            refCell.lba() = vector(0, 0, 0);
            refCell.rba() = vector(cellSizes[0]/double(Nx), 0, 0);
            refCell.lta() = vector(0, cellSizes[1]/double(Ny), 0);
            refCell.rta() = vector(cellSizes[0]/double(Nx), cellSizes[1]/double(Ny), 0);
            refCell.lbf() = vector(0, 0, cellSizes[2]/double(Nz));
            refCell.rbf() = vector(cellSizes[0]/double(Nx), 0, cellSizes[2]/double(Nz));
            refCell.ltf() = vector(0, cellSizes[1]/double(Ny), cellSizes[2]/double(Nz));
            refCell.rtf() = vector(cellSizes[0]/double(Nx), cellSizes[1]/double(Ny), cellSizes[2]/double(Nz));

            forAllCells(alpha, i, j, k)
            {

                scalar exactVol = 0;

                for (int aux1 = 0; aux1 < Nx; aux1++)
                {
                    for (int aux2 = 0; aux2 < Ny; aux2++)
                    {
                        for (int aux3 = 0; aux3 < Nz; aux3++)
                        {
                            vertexVector cell = vertex(i,j,k).lba()
                                + (double(aux1)) * vector(cellSizes[0]/double(Nx), 0, 0)
                                + (double(aux2)) * vector(0, cellSizes[1]/double(Ny), 0)
                                + (double(aux3)) * vector(0, 0, cellSizes[2]/double(Nz))
                                + refCell;

                            for (int l = 0; l < 8; l++)
                            {
                                vector v = cell[l];
                                exactVol += (1/(double(Nx * Ny * Nz) * cellSizes[0] * cellSizes[1] * cellSizes[2])) * 0.125
                                        * Foam::min(cellSizes[0] * cellSizes[1] * cellSizes[2], Foam::max(0,
                                        10 * (Foam::sqr(0.15) - Foam::sqr(Foam::mag(v - vector(0.35,0.35,0.35))))));
                            }
                        }
                    }
                }

                exactVol = Foam::min(1, exactVol);
                exactVol = Foam::max(0, exactVol);

                L1Error += cv(i,j,k) * Foam::mag(exactVol - alpha(i,j,k));
                Volume += cv(i,j,k) * alpha(i,j,k);
                LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);

                if (LocalBoundError > BoundError)
                    BoundError = LocalBoundError;
            }
        }
        else
        {
            forAllCells(alpha, i, j, k)
            {
                scalar exactVol = 0;

                for (int l = 0; l < 8; l++)
                {
                    vector v = vertex(i,j,k)[l];
                    exactVol += 0.125 * (Foam::mag(v - vector(0.35,0.35,0.35)) < 0.15);
                }

                L1Error += cv(i,j,k) * Foam::mag(exactVol - alpha(i,j,k));
                Volume += cv(i,j,k) * alpha(i,j,k);
                LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);

            if (LocalBoundError > BoundError)
                BoundError = LocalBoundError;
            }
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
    }

    return true;
}

}

}

}

}
