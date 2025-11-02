#include "curvatureDisk.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

using ::Foam::sqr;
using ::Foam::sqrt;
using ::Foam::min;
using ::Foam::max;
using ::Foam::asin;

defineTypeNameAndDebug(curvatureDisk, 0);

addToRunTimeSelectionTable
(
    functionObject,
    curvatureDisk,
    dictionary
);

bool curvatureDisk::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        const fvMesh& fvMsh = alpha.fvMsh();

        const meshField<vertexVector,colocated>& vertex =
            fvMsh.metrics<colocated>().vertexCenters();

        alpha = Zero;

        forAllCells(alpha, i, j, k)
        {
            const labelVector ijk(i,j,k);

            vector Min = min(vertex(ijk).lba(), vertex(ijk).rtf());
            vector Max = max(vertex(ijk).lba(), vertex(ijk).rtf());

            scalar x0 = Min.x();
            scalar x1 = Max.x();
            scalar y0 = Min.y();
            scalar y1 = Max.y();

            // Flip

            if (x0 < -1e-12)
            {
                scalar t = x0;
                x0 = -x1;
                x1 = -t;
            }

            if (y0 < -1e-12)
            {
                scalar t = y0;
                y0 = -y1;
                y1 = -t;
            }

            const scalar totalVol = (x1 - x0)*(y1 - y0);

            if (sqr(x1) + sqr(y1) < sqr(R_))
            {
                alpha(ijk) = 1.0;
            }
            else if (sqr(x0) + sqr(y1) < sqr(R_))
            {
                if (sqr(x1) + sqr(y0) < sqr(R_))
                {
                    const scalar aux = sqrt(sqr(R_) - sqr(y1));

                    const scalar vol =
                        0.5
                      * (
                            sqr(R_)*asin(x1/R_)
                          - sqr(R_)*asin(aux/R_)
                          + x1*sqrt(sqr(R_) - sqr(x1))
                          - aux*sqrt(sqr(R_) - sqr(aux))
                        )
                      - y0*(x1 - aux)
                      + (y1 - y0)*(aux - x0);

                    alpha(ijk) = vol/totalVol;
                }
                else
                {
                    const scalar vol =
                        0.5
                      * (
                            sqr(R_)*asin(y1/R_)
                          - sqr(R_)*asin(y0/R_)
                          + y1*sqrt(sqr(R_) - sqr(y1))
                          - y0*sqrt(sqr(R_) - sqr(y0))
                        )
                      - x0*(y1 - y0);

                    alpha(ijk) = vol/totalVol;
                }
            }
            else if (sqr(x1) + sqr(y0) < sqr(R_))
            {
                const scalar vol =
                    0.5
                  * (
                        sqr(R_)*asin(x1/R_)
                      - sqr(R_)*asin(x0/R_)
                      + x1*sqrt(sqr(R_) - sqr(x1))
                      - x0*sqrt(sqr(R_) - sqr(x0))
                    )
                  - y0*(x1 - x0);

                alpha(ijk) = vol/totalVol;
            }
            else if (sqr(x0) + sqr(y0) < sqr(R_))
            {
                const scalar aux = sqrt(sqr(R_) - sqr(y0));

                const scalar vol =
                    0.5
                  * (
                        sqr(R_)*asin(aux/R_)
                      - sqr(R_)*asin(x0/R_)
                      + aux*sqrt(sqr(R_) - sqr(aux))
                      - x0*sqrt(sqr(R_) - sqr(x0))
                    )
                  - y0*(aux - x0);

                alpha(ijk) = vol/totalVol;
            }
        }

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool curvatureDisk::end()
{
    colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    colocatedScalarField& kappa =
        runTime_.lookupObjectRef<colocatedScalarField>("kappa");

    scalar error = 0;
    label count = 0;

    forAllCells(kappa, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if
        (
            alpha(ijk) >       vof::threshold
         && alpha(ijk) < 1.0 - vof::threshold
        )
        {
            error += sqr(kappa(ijk) - 1/R_);
            count++;
        }
    }

    reduce(error, sumOp<scalar>());
    reduce(count, sumOp<label>());

    error = R_*sqrt(error/scalar(count));

    Info<< "Real curvature: " << 1/R_ << endl
        << "Curvature error: " << error << endl;

    return true;
}

}

}

}

}
