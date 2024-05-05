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
    cAlpha = Zero;

    forAllCells(alpha_, i, j, k)
        for (int d = 0; d < 3; d++)
            cAlpha(i,j,k) +=
                0.5 /3.0*alpha_(i,j,k)
              + 0.25/3.0*alpha_(lowerNei(i,j,k,d))
              + 0.25/3.0*alpha_(upperNei(i,j,k,d));

    cAlpha.correctBoundaryConditions();

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

    normal.correctBoundaryConditions();

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
