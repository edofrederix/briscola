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
    normalScheme(vf, dict)
{}

Youngs::Youngs(const Youngs& s)
:
    normalScheme(s)
{}

Youngs::~Youngs()
{}

tmp<colocatedVectorField> Youngs::operator()()
{
    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(ex::grad(vf_.alpha()))
    );

    const colocatedScalarDirection& alpha = vf_.alpha()[0][0];
    colocatedVectorDirection& n = tn.ref()[0][0];

    forAllCells(n, i, j, k)
    {
        if (alpha(i,j,k) > 1e-6 && alpha(i,j,k) < 1-1e-6)
        {
            n(i,j,k) /= Foam::mag(n(i,j,k));
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
