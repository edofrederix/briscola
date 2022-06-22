#include "geometricGrading.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace gradings
{

defineTypeNameAndDebug(geometricGrading, 0);
addToRunTimeSelectionTable(grading, geometricGrading, dictionary);

geometricGrading::geometricGrading(const brick& b)
:
    grading(b)
{
    Istream& is = b.dict().lookup("grading");

    word dummy;

    is >> dummy;
    is >> f_;

    for (label cmpt = 0; cmpt < 3; cmpt++)
        if (f_[cmpt] <= 0.0)
            FatalErrorInFunction
                << "Component " << cmpt << " of geometric grading in brick "
                << b.num() << " is invalid." << endl
                << abort(FatalError);
}

tmp<scalarField> geometricGrading::operator()
(
    const scalarField& f,
    const label e
) const
{
    const label cmpt = e/4;

    // Number of cells

    const scalar n = scalar(this->b_.N()[cmpt]);

    if (n == 1.0)
    {
        return tmp<scalarField>(new scalarField(f));
    }
    else
    {
        const scalar g = Foam::pow(f_[cmpt], 1.0/(n-1.0));

        if (g == 1.0)
        {
            return tmp<scalarField>(new scalarField(f));
        }
        else
        {
            const scalar d = (g-1)/(g*(Foam::pow(g,n-1)-1));

            return d*g*(Foam::pow(g,f*(n-1))-1)/(g-1);
        }
    }
}

}

}

}
