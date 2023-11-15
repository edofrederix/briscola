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
    briscolaFunctionObject(name, runTime, dict),
    T_(3.0),
    TotalVolume_(0.0),
    BoundError_(0.0),
    Nx_(2),
    Ny_(2),
    Nz_(2)
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

        const mesh& mh = fvMsh.msh();

        forAllCells(alpha, i, j, k)
        {
            vector ccxy = cc(i,j,k);

            U(i,j,k) = vector
                (
                  2 * Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.x()))
                    * Foam::sin(2 * pi * ccxy.y())
                    * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.y()))
                    * Foam::sin(2 * pi * ccxy.x())
                    * Foam::sin(2 * pi * ccxy.z()),
                  - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.z()))
                    * Foam::sin(2 * pi * ccxy.y())
                    * Foam::sin(2 * pi * ccxy.x())
                );
        }

        if (mh.castable<uniformMesh>())
        {
            const uniformMesh& uniMh = mh.cast<uniformMesh>();
            vector cellSizes = uniMh.cellSize();

            vertexVector refCell;

            refCell.lba() = vector(0, 0, 0);
            refCell.rba() = vector(cellSizes[0]/double(Nx_), 0, 0);
            refCell.lta() = vector(0, cellSizes[1]/double(Ny_), 0);
            refCell.rta() = vector(cellSizes[0]/double(Nx_), cellSizes[1]/double(Ny_), 0);
            refCell.lbf() = vector(0, 0, cellSizes[2]/double(Nz_));
            refCell.rbf() = vector(cellSizes[0]/double(Nx_), 0, cellSizes[2]/double(Nz_));
            refCell.ltf() = vector(0, cellSizes[1]/double(Ny_), cellSizes[2]/double(Nz_));
            refCell.rtf() = vector(cellSizes[0]/double(Nx_), cellSizes[1]/double(Ny_), cellSizes[2]/double(Nz_));

            forAllCells(alpha, i, j, k)
            {

                alpha(i,j,k) = 0;

                for (int aux1 = 0; aux1 < Nx_; aux1++)
                {
                    for (int aux2 = 0; aux2 < Ny_; aux2++)
                    {
                        for (int aux3 = 0; aux3 < Nz_; aux3++)
                        {
                            vertexVector cell = vertex(i,j,k).lba()
                                + (double(aux1)) * vector(cellSizes[0]/double(Nx_), 0, 0)
                                + (double(aux2)) * vector(0, cellSizes[1]/double(Ny_), 0)
                                + (double(aux3)) * vector(0, 0, cellSizes[2]/double(Nz_))
                                + refCell;

                            for (int l = 0; l < 8; l++)
                            {
                                vector v = cell[l];
                                alpha(i,j,k) += (1/(double(Nx_ * Ny_ * Nz_) * cellSizes[0] * cellSizes[1] * cellSizes[2])) * 0.125
                                        * Foam::min(cellSizes[0] * cellSizes[1] * cellSizes[2], Foam::max(0,
                                        10 * (Foam::sqr(0.15) - Foam::sqr(Foam::mag(v - vector(0.35,0.35,0.35))))));
                            }
                        }
                    }
                }

                alpha(i,j,k) = Foam::min(1, alpha(i,j,k));
                alpha(i,j,k) = Foam::max(0, alpha(i,j,k));

                TotalVolume_ += cv(i,j,k) * alpha(i,j,k);
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
                TotalVolume_ += cv(i,j,k) * alpha(i,j,k);
            }
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
    const scalar pi = constant::mathematical::pi;
    scalar LocalBoundError;

    colocatedVectorField& U =
        runTime_.lookupObjectRef<colocatedVectorField>("U");

    colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const fvMesh& fvMsh = U.fvMsh();

    const colocatedVectorField& cc =
        fvMsh.metrics<colocated>().cellCenters();

    const colocatedScalarField& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    forAllCells(U, i, j, k)
    {
        vector ccxy = cc(i,j,k);

        U(i,j,k) = vector
            (
                2 * Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value())/ T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.x()))
                    * Foam::sin(2 * pi * ccxy.y())
                    * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.y()))
                    * Foam::sin(2 * pi * ccxy.x())
                    * Foam::sin(2 * pi * ccxy.z()),
                - Foam::cos(pi * (runTime_.time().value() + 0.5 * runTime_.deltaT().value()) / T_)
                    * Foam::sqr(Foam::sin(pi * ccxy.z()))
                    * Foam::sin(2 * pi * ccxy.y())
                    * Foam::sin(2 * pi * ccxy.x())
            );

        LocalBoundError = cv(i,j,k) * Foam::max(-alpha(i,j,k), alpha(i,j,k)-1);
        if (LocalBoundError > BoundError_)
            BoundError_ = LocalBoundError;

    }

    reduce(BoundError_, maxOp<scalar>());

    colocatedFaceScalarField& phi =
        runTime_.lookupObjectRef<colocatedFaceScalarField>("phi");

    U.correctBoundaryConditions();

    phi = ex::faceFlux(U);

    return true;
}

bool initialCondition::end()
{
    if (Foam::mag(runTime_.time().value() - T_) < 1e-12)
    {

        scalar L1Error = 0;
        scalar Volume = 0;
        scalar LocalBoundError;

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const fvMesh& fvMsh = alpha.fvMsh();

        const colocatedScalarField& cv =
            fvMsh.metrics<colocated>().cellVolumes();

        const meshField<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters();

        const mesh& mh = fvMsh.msh();

        if (mh.castable<uniformMesh>())
        {
            const uniformMesh& uniMh = mh.cast<uniformMesh>();
            vector cellSizes = uniMh.cellSize();

            vertexVector refCell;

            refCell.lba() = vector(0, 0, 0);
            refCell.rba() = vector(cellSizes[0]/double(Nx_), 0, 0);
            refCell.lta() = vector(0, cellSizes[1]/double(Ny_), 0);
            refCell.rta() = vector(cellSizes[0]/double(Nx_), cellSizes[1]/double(Ny_), 0);
            refCell.lbf() = vector(0, 0, cellSizes[2]/double(Nz_));
            refCell.rbf() = vector(cellSizes[0]/double(Nx_), 0, cellSizes[2]/double(Nz_));
            refCell.ltf() = vector(0, cellSizes[1]/double(Ny_), cellSizes[2]/double(Nz_));
            refCell.rtf() = vector(cellSizes[0]/double(Nx_), cellSizes[1]/double(Ny_), cellSizes[2]/double(Nz_));

            forAllCells(alpha, i, j, k)
            {

                scalar exactVol = 0;

                for (int aux1 = 0; aux1 < Nx_; aux1++)
                {
                    for (int aux2 = 0; aux2 < Ny_; aux2++)
                    {
                        for (int aux3 = 0; aux3 < Nz_; aux3++)
                        {
                            vertexVector cell = vertex(i,j,k).lba()
                                + (double(aux1)) * vector(cellSizes[0]/double(Nx_), 0, 0)
                                + (double(aux2)) * vector(0, cellSizes[1]/double(Ny_), 0)
                                + (double(aux3)) * vector(0, 0, cellSizes[2]/double(Nz_))
                                + refCell;

                            for (int l = 0; l < 8; l++)
                            {
                                vector v = cell[l];
                                exactVol += (1/(double(Nx_ * Ny_* Nz_) * cellSizes[0] * cellSizes[1] * cellSizes[2])) * 0.125
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

                if (LocalBoundError > BoundError_)
                    BoundError_ = LocalBoundError;
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

            if (LocalBoundError > BoundError_)
                BoundError_ = LocalBoundError;
            }
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
    }

    return true;
}

}

}

}

}
