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
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

Youngs::Youngs(const Youngs& s)
:
    normalScheme(s),
    threshold_(s.threshold_)
{}

Youngs::~Youngs()
{}

tmp<colocatedVectorField> Youngs::operator()()
{
    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(ex::grad(vf_.alpha()))
    );

    colocatedVectorDirection& n = tn.ref()[0][0];

    forAllCells(n, i, j, k)
    {
        const scalar S = Foam::mag(n(i,j,k));

        if (S > threshold_)
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
