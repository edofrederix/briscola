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
    gradThreshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

Youngs::Youngs(const Youngs& s)
:
    normalScheme(s),
    gradThreshold_(s.gradThreshold_)
{}

Youngs::~Youngs()
{}

void Youngs::correct()
{
    colocatedVectorField& n = *this;

    n = ex::grad(alpha_);

    forAllCells(n, i, j, k)
    {
        const scalar S = Foam::mag(n(i,j,k));

        if (S > gradThreshold_)
        {
            n(i,j,k) /= S;
        }
        else
        {
            n(i,j,k) = Zero;
        }
    }

    n[0].correctBoundaryConditions();
}

}

}

}
