#include "CV.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceTensionScheme.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(CV, 0);
addToRunTimeSelectionTable(curvatureScheme, CV, dictionary);

CV::CV
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    curvatureScheme(fvMsh, dict, normal, alpha)
{}

CV::CV(const CV& s)
:
    curvatureScheme(s)
{}

CV::~CV()
{}

void CV::correct()
{
    colocatedScalarField& kappa = *this;

    colocatedScalarField cAlpha("cAlpha", fvMsh_);

    scalar kernel[3][3] =
    {
        {1.0/16.0, 2.0/16.0, 1.0/16.0},
        {2.0/16.0, 4.0/16.0, 2.0/16.0},
        {1.0/16.0, 2.0/16.0, 1.0/16.0}
    };

    cAlpha = Zero;

    forAllCells(alpha_, i, j, k)
    {
        for (int aux1 = -1; aux1 <= 1; aux1++)
        {
            for (int aux2 = -1; aux2 <= 1; aux2++)
            {
                cAlpha(i,j,k) +=
                    kernel[aux1 + 1][aux2 + 1]
                  * alpha_(i+aux1,j+aux2,k);
            }
        }
    }

    colocatedVectorField normal(ex::grad(cAlpha));

    forAllCells(alpha_, i, j, k)
    {
        scalar S = Foam::mag(normal(i,j,k));

        if (S > 1e-8)
        {
            normal(i,j,k) /= S;
        }
        else
        {
            normal(i,j,k) = Zero;
        }
    }

    const colocatedScalarField kappa2(ex::div(ex::faceFlux(normal)()));

    kappa = Zero;

    forAllCells(alpha_, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) >     vof::threshold)
         && (alpha_(i,j,k) < 1 - vof::threshold)
        )
        {
            int count = 0;

            for (int aux1 = -1; aux1 <= 1; aux1++)
            {
                for (int aux2 = -1; aux2 <= 1; aux2++)
                {
                    if
                    (
                        (alpha_(i+aux1,j+aux2,k) >     vof::threshold)
                     && (alpha_(i+aux1,j+aux2,k) < 1 - vof::threshold)
                    )
                    {
                        count++;
                        kappa(i,j,k) += kappa2(i+aux1,j+aux2,k);
                    }
                }
            }

            kappa(i,j,k) /= scalar(count);
        }
    }
}

}

}

}
