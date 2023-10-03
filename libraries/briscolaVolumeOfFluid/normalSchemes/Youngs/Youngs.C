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

Youngs::Youngs(const vof& vf, const dictionary& dict)
:
    normalScheme(vf, dict),
    gradThreshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

Youngs::Youngs(const Youngs& s)
:
    normalScheme(s),
    gradThreshold_(s.gradThreshold_)
{}

Youngs::~Youngs()
{}

tmp<colocatedVectorField> Youngs::operator()()
{
    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(ex::grad(vf_.alpha()))
    );

    colocatedVectorField& n = tn.ref();

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

    return tn;
}

}

}

}
