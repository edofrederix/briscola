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
    curvatureScheme(fvMsh, dict, normal, alpha),
    YoungsNormal_(fvMsh, alpha)
{}

CV::CV(const CV& s)
:
    curvatureScheme(s),
    YoungsNormal_(s.YoungsNormal_)
{}

CV::~CV()
{}

void CV::correct()
{
    colocatedScalarField& kappa = *this;

    YoungsNormal_.correct();

    // Curvature

    kappa = - ex::div(ex::faceFlux(YoungsNormal_)());

    // Clip for non-interfacial cells

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
