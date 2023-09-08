#include "LSGIR.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(LSGIR, 0);
addToRunTimeSelectionTable(normalScheme, LSGIR, dictionary);

LSGIR::LSGIR(const vof& vf, const dictionary& dict)
:
    normalScheme(vf, dict),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

LSGIR::LSGIR(const LSGIR& s)
:
    normalScheme(s),
    threshold_(s.threshold_)
{}

LSGIR::~LSGIR()
{}

tmp<colocatedVectorField> LSGIR::operator()()
{

    colocatedScalarField aux = vf_.alpha();
    colocatedScalarDirection alpha = aux[0][0];

    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(
            "normal",
            vf_.fvMsh()
        )
    );

    colocatedVectorDirection& n = tn.ref()[0][0];

    const meshDirection<vector,colocated>& centers =
            vf_.fvMsh().template
            metrics<colocated>().cellCenters()[0][0];


    double Aaux[26][3];
    double baux[26];
    int index;
    scalarList weights(26);


    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > 1e-12) && (alpha(i,j,k) < 1 - 1e-12))
        {
            tensor A = Zero;
            vector b = Zero;

            for (int aux1 = 0; aux1 < 3; aux1++)
            {
                for (int aux2 = 0; aux2 < 3; aux2++)
                {
                    for (int aux3 = 0; aux3 < 3; aux3++)
                    {
                        if ((aux1 != 1) || (aux2 != 1) || (aux3 != 1))
                        {
                            index = aux1 + 3 * aux2 + 9 * aux3;
                            if (index > 13)
                                index--;
                            weights[index] = (1.0) / Foam::sqr(Foam::mag(centers(i,j,k) - centers(i+aux1-1,j+aux2-1,k+aux3-1)));
                            baux[index] = weights[index] * (alpha(i+aux1-1,j+aux2-1,k+aux3-1) - alpha(i,j,k));
                            Aaux[index][0] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[0] - centers(i,j,k)[0]);
                            Aaux[index][1] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[1] - centers(i,j,k)[1]);
                            Aaux[index][2] = weights[index] * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[2] - centers(i,j,k)[2]);
                        }
                    }
                }
            }

            for (int aux3 = 0; aux3 < 26; aux3++)
            {
                for (int aux1 = 0; aux1 < 3; aux1++)
                {
                    for (int aux2 = 0; aux2 < 3; aux2++)
                    {
                        A[3 * aux1 + aux2] += Aaux[aux3][aux1] * Aaux[aux3][aux2];
                    }
                    b[aux1] += Aaux[aux3][aux1] * baux[aux3];
                }
            }

            n(i,j,k) = A.inv() & b;
            const scalar S = Foam::mag(n(i,j,k));
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
