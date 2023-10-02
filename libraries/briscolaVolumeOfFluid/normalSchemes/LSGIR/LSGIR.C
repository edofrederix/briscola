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
    threshold_(vf.threshold())
{
    createBoundaryType();
}

LSGIR::LSGIR(const LSGIR& s)
:
    normalScheme(s),
    threshold_(s.threshold_)
{
    createBoundaryType();
}

LSGIR::~LSGIR()
{}

void LSGIR::createBoundaryType()
{
    const faceLabel faceType =
        vf_.fvMsh().msh().patchType();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                bool aux = true;

                if ((vf_.fvMsh()[0].l() > 1) && (i != 1))
                {
                    aux = (i == 0) ? (aux && (faceType.left() > 0)) : (aux && (faceType.right() > 0));
                }

                if ((vf_.fvMsh()[0].m() > 1) && (j != 1))
                {
                    aux = (j == 0) ? (aux && (faceType.bottom() > 0)) : (aux && (faceType.top() > 0));
                }

                if ((vf_.fvMsh()[0].n() > 1) && (k != 1))
                {
                    aux = (k == 0) ? (aux && (faceType.aft() > 0)) : (aux && (faceType.fore() > 0));
                }

                boundaryType_[i][j][k] = aux;

            }
        }
    }
}

tmp<colocatedVectorField> LSGIR::operator()()
{

    /*

    +++++++++++++++++++++++++++++++++++++++++++++++

    Reconstruc the interface normal using the version
    of LSGIR presented in Lopez (2022).

    +++++++++++++++++++++++++++++++++++++++++++++++

    */

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

    int index;
    scalar weight;

    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > threshold_) && (alpha(i,j,k) < 1 - threshold_))
        {
            double Aaux[26][3] = {0};
            double baux[26] = {0};
            tensor A = Zero;
            vector b = Zero;

            for (int aux1 = 0; aux1 < 3; aux1++)
            {
                for (int aux2 = 0; aux2 < 3; aux2++)
                {
                    for (int aux3 = 0; aux3 < 3; aux3++)
                    {
                        bool interiorNode = (((i+aux1-1) >= n.I().left()) &&
                            ((i+aux1-1) < n.I().right()) &&
                            ((j+aux2-1) >= n.I().bottom()) &&
                            ((j+aux2-1) < n.I().top()) &&
                            ((k+aux3-1) >= n.I().aft()) &&
                            ((k+aux3-1) < n.I().fore()));

                        if (((aux1 != 1) || (aux2 != 1) || (aux3 != 1)) &&
                            (interiorNode || boundaryType_[aux1][aux2][aux3]))
                            {
                                index = aux1 + 3 * aux2 + 9 * aux3;
                                if (index > 13)
                                    index--;
                                weight = (1.0) / Foam::pow(Foam::mag(centers(i,j,k) - centers(i+aux1-1,j+aux2-1,k+aux3-1)),1.5);
                                baux[index] = weight * (alpha(i+aux1-1,j+aux2-1,k+aux3-1) - alpha(i,j,k));
                                Aaux[index][0] = weight * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[0] - centers(i,j,k)[0]);
                                Aaux[index][1] = weight * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[1] - centers(i,j,k)[1]);
                                Aaux[index][2] = weight * (centers(i+aux1-1,j+aux2-1,k+aux3-1)[2] - centers(i,j,k)[2]);
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
