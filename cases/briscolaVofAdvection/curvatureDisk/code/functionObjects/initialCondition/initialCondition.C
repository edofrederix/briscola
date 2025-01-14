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
    R_(0.5)
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

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const fvMesh& fvMsh = U.fvMsh();

        const meshField<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters();

        forAllCells(alpha, i, j, k)
        {
            U(i,j,k) = vector
            (
                0,
                0,
                0
            );
        }

        forAllCells(alpha, i, j, k)
        {
            scalar x0 = vertex(i,j,k).lba()[0];
            scalar x1 = vertex(i,j,k).rba()[0];
            scalar y0 = vertex(i,j,k).lba()[1];
            scalar y1 = vertex(i,j,k).lta()[1];

            if (x0 < -1e-12)
            {
                scalar aux = x0;
                x0 = -x1;
                x1 = -aux;
            }

            if (y0 < -1e-12)
            {
                scalar aux = y0;
                y0 = -y1;
                y1 = -aux;
            }

            scalar totalVol = (x1 - x0) * (y1 - y0);

            if (Foam::sqr(x1) + Foam::sqr(y1) < Foam::sqr(R_))
            {
                alpha(i,j,k) = 1;
            }
            else if (Foam::sqr(x0) + Foam::sqr(y1) < Foam::sqr(R_))
            {
                if (Foam::sqr(x1) + Foam::sqr(y0) < Foam::sqr(R_))
                {
                    scalar aux = Foam::sqrt(Foam::sqr(R_) - Foam::sqr(y1));

                    scalar vol = 0.5 *
                        (
                            Foam::sqr(R_) * Foam::asin(x1 / R_)
                          - Foam::sqr(R_) * Foam::asin(aux / R_)
                          + x1 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(x1))
                          - aux * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(aux))
                        )
                      - y0 * (x1 - aux)
                      + (y1 - y0) * (aux - x0);

                    alpha(i,j,k) = vol / totalVol;

                }
                else
                {
                    scalar vol = 0.5 *
                        (
                            Foam::sqr(R_) * Foam::asin(y1 / R_)
                          - Foam::sqr(R_) * Foam::asin(y0 / R_)
                          + y1 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(y1))
                          - y0 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(y0))
                        )
                      - x0 * (y1 - y0);

                    alpha(i,j,k) = vol / totalVol;
                }
            }
            else if (Foam::sqr(x1) + Foam::sqr(y0) < Foam::sqr(R_))
            {
                scalar vol = 0.5 *
                    (
                        Foam::sqr(R_) * Foam::asin(x1 / R_)
                      - Foam::sqr(R_) * Foam::asin(x0 / R_)
                      + x1 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(x1))
                      - x0 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(x0))
                    )
                  - y0 * (x1 - x0);

                alpha(i,j,k) = vol / totalVol;
            }
            else if (Foam::sqr(x0) + Foam::sqr(y0) < Foam::sqr(R_))
            {
                scalar aux = Foam::sqrt(Foam::sqr(R_) - Foam::sqr(y0));

                scalar vol = 0.5 *
                    (
                        Foam::sqr(R_) * Foam::asin(aux / R_)
                      - Foam::sqr(R_) * Foam::asin(x0 / R_)
                      + aux * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(aux))
                      - x0 * Foam::sqrt(Foam::sqr(R_) - Foam::sqr(x0))
                    )
                  - y0 * (aux - x0);

                alpha(i,j,k) = vol / totalVol;

            }
            else
            {
                alpha(i,j,k) = 0;
            }

        }

        alpha.correctBoundaryConditions();

        colocatedFaceScalarField& phi =
            runTime_.lookupObjectRef<colocatedFaceScalarField>("phi");

        U.correctBoundaryConditions();

        phi = ex::faceFlux(U);
    }

    return true;
}

bool initialCondition::end()
{
    scalar error = 0;
    int count = 0;

    colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    colocatedScalarField& kappa =
        runTime_.lookupObjectRef<colocatedScalarField>("kappa");

    forAllCells(kappa, i, j, k)
    {

        if
        (
            alpha(i,j,k) > 1e-8
            && alpha(i,j,k) < (1.0 - 1e-8)
        )
        {
            error += Foam::sqr(kappa(i,j,k) - 1/R_);
            count++;
        }
    }

    reduce(error, sumOp<scalar>());
    reduce(count, sumOp<int>());
    error = R_ * Foam::sqrt(error/double(count));
    Info << "Real curvature: " << 1/R_ << endl
         << "Curvature error: " << error << endl;

    return true;
}

}

}

}

}
