#include "pipeFlowIBM.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "meshFields.H"

namespace Foam
{

using constant::mathematical::pi;

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(pipeFlowIBM, 0);

addToRunTimeSelectionTable
(
    functionObject,
    pipeFlowIBM,
    dictionary
);

bool pipeFlowIBM::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        if (runTime_.foundObject<colocatedVectorField>("U"))
        {
            colocatedVectorField& U =
                runTime_.lookupObjectRef<colocatedVectorField>("U");

            const colocatedVectorField& cc =
                U.fvMsh().metrics<colocated>().cellCenters();

            forAllCells(U, i, j, k)
            {
                const scalar x = cc(i,j,k).x();
                const scalar y = cc(i,j,k).y();
                const scalar z = cc(i,j,k).z();

                const scalar r = Foam::sqrt(Foam::sqr(y)+Foam::sqr(z));
                const scalar theta = Foam::atan2(y,z);

                U(i,j,k).x() =
                    5300.0/360.0
                  * (3.0*(1.0-Foam::sqr(r)) + Foam::sin(r*8.0*pi)/4.0);

                U(i,j,k).y() =
                    5300.0/360.0*Foam::sin(r*pi*2.0)
                  * Foam::sin(x*2.0)*Foam::cos(theta)/4.0;

                U(i,j,k).z() =
                    5300.0/360.0*Foam::sin(r*pi*2.0)
                  * Foam::sin(x*2.0)*Foam::sin(theta)/4.0;
            }

            U.correctBoundaryConditions();
        }
        else
        {
            staggeredScalarField& U =
                runTime_.lookupObjectRef<staggeredScalarField>("U");

            const fvMesh& fvMsh = U.fvMsh();

            const tensor base =
                fvMsh.msh().cast<rectilinearMesh>().base();

            const staggeredVectorField& cc =
                fvMsh.metrics<staggered>().cellCenters();

            forAllCells(U, d, i, j, k)
            {
                const scalar x = cc(d,i,j,k).x();
                const scalar y = cc(d,i,j,k).y();
                const scalar z = cc(d,i,j,k).z();

                const scalar r = Foam::sqrt(Foam::sqr(y)+Foam::sqr(z));
                const scalar theta = Foam::atan2(y,z);

                const vector u
                (
                    5300.0/360.0
                  * (3.0*(1.0-Foam::sqr(r)) + Foam::sin(r*8.0*pi)/4.0),
                    5300.0/360.0*Foam::sin(r*pi*2.0)
                  * Foam::sin(x*2.0)*Foam::cos(theta)/4.0,
                    5300.0/360.0*Foam::sin(r*pi*2.0)
                  * Foam::sin(x*2.0)*Foam::sin(theta)/4.0
                );

                U(d,i,j,k) = (base & u)[d];
            }

            U.correctBoundaryConditions();
        }
    }

    return true;
}

}

}

}

}