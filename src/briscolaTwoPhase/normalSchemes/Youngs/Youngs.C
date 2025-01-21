#include "Youngs.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(Youngs, 0);
addToRunTimeSelectionTable(normalScheme, Youngs, dictionary);

Youngs::Youngs
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dict, alpha),
    nSmooth_(dict.lookupOrDefault<scalar>("nSmooth", 1))
{}

Youngs::Youngs
(
    const fvMesh& fvMsh,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dictionary::null, alpha),
    nSmooth_(1)
{}

Youngs::Youngs(const Youngs& s)
:
    normalScheme(s),
    nSmooth_(s.nSmooth_)
{}

Youngs::~Youngs()
{}

void Youngs::correct()
{
    colocatedVectorField& n = *this;

    colocatedScalarField alpha(alpha_);

    for (int i = 0; i < nSmooth_; i++)
    {
        // Jacobi sweep (Gauss-Seidel will make every cell interfacial)

        colocatedScalarField alpha0(alpha);

        forAllCells(alpha, i, j, k)
            for (int d = 0; d < 3; d++)
                alpha(i,j,k) =
                    0.5 *alpha0(i,j,k)
                  + 0.25*alpha0(lowerNeighbor(i,j,k,d))
                  + 0.25*alpha0(upperNeighbor(i,j,k,d));
    }

    const colocatedFaceScalarField& fa =
        fvMsh_.metrics<colocated>().faceAreas();

    n = ex::reconstruct(ex::faceGrad(alpha)*fa);

    forAllCells(n, i, j, k)
        if (Foam::mag(n(i,j,k)))
            n(i,j,k) /= Foam::mag(n(i,j,k));

    n.correctBoundaryConditions();
}

}

}

}
