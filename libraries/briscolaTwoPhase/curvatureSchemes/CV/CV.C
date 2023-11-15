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

    colocatedScalarField cAlpha("alpha", fvMsh_);

    scalar kernel[3] =
    {
        0.25, 0.5, 0.25
    };

    cAlpha = Zero;

    forAllCells(alpha_, i, j, k)
    {
        for (int aux1 = -1; aux1 <= 1; aux1++)
        {
            for (int aux2 = -1; aux2 <= 1; aux2++)
            {
                for (int aux3 = -1; aux3 <= 1; aux3++)
                {
                    cAlpha(i,j,k) +=
                        kernel[aux1 + 1]
                      * kernel[aux2 + 1]
                      * kernel[aux3 + 1]
                      * alpha_(i+aux1,j+aux2,k+aux3);
                }
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

    kappa = - ex::div(ex::faceFlux(normal)());

    forAllCells(alpha_, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) <=     vof::threshold)
         || (alpha_(i,j,k) >= 1 - vof::threshold)
        )
        {
            kappa(i,j,k) = 0;
        }
    }

    kappa.correctBoundaryConditions();
}

}

}

}
