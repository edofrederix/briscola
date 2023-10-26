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

void LSGIR::createBoundaryType()
{
    const faceLabel& faceType =
        vf_.fvMsh().msh().facePatchType();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                bool aux = true;

                if ((vf_.fvMsh()[0].l() > 1) && (i != 1))
                {
                    aux =
                        (i == 0)
                      ? (aux && (faceType.left()  > 0))
                      : (aux && (faceType.right() > 0));
                }

                if ((vf_.fvMsh()[0].m() > 1) && (j != 1))
                {
                    aux =
                        (j == 0)
                      ? (aux && (faceType.bottom() > 0))
                      : (aux && (faceType.top()    > 0));
                }

                if ((vf_.fvMsh()[0].n() > 1) && (k != 1))
                {
                    aux =
                        (k == 0)
                      ? (aux && (faceType.aft()  > 0))
                      : (aux && (faceType.fore() > 0));
                }

                boundaryType_[i][j][k] = aux;
            }
        }
    }
}

LSGIR::LSGIR(const vof& vf, const dictionary& dict)
:
    normalScheme(vf, dict)
{
    createBoundaryType();
}

LSGIR::LSGIR(const LSGIR& s)
:
    normalScheme(s)
{
    createBoundaryType();
}

LSGIR::~LSGIR()
{}

tmp<colocatedVectorField> LSGIR::operator()()
{
    // Reconstruct the interface normal using the version of LSGIR presented in
    // Lopez (2022)

    const colocatedScalarField& alpha = vf_.alpha();

    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField
        (
            "normal",
            vf_.fvMsh()
        )
    );

    colocatedVectorField& n = tn.ref();

    const meshField<vector,colocated>& centers =
        vf_.fvMsh().template metrics<colocated>().cellCenters();

    int index;
    scalar weight;

    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > vof::threshold) && (alpha(i,j,k) < 1 - vof::threshold))
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
                        bool interiorNode =
                            (
                                (i+aux1-1) >= n.I().left()
                             && (i+aux1-1) <  n.I().right()
                             && (j+aux2-1) >= n.I().bottom()
                             && (j+aux2-1) <  n.I().top()
                             && (k+aux3-1) >= n.I().aft()
                             && (k+aux3-1) <  n.I().fore()
                            );

                        if
                        (
                            (aux1 != 1 || aux2 != 1 || aux3 != 1)
                         && (interiorNode || boundaryType_[aux1][aux2][aux3])
                        )
                        {
                            index = aux1 + 3 * aux2 + 9 * aux3;

                            if (index > 13)
                                index--;

                            weight =
                                1.0/Foam::pow
                                (
                                    Foam::mag
                                    (
                                        centers(i,j,k)
                                      - centers(i+aux1-1,j+aux2-1,k+aux3-1)
                                    ),
                                    1.5
                                );

                            baux[index] =
                                weight *
                                (
                                    alpha(i+aux1-1,j+aux2-1,k+aux3-1)
                                  - alpha(i,j,k)
                                );

                            Aaux[index][0] =
                                weight *
                                (
                                    centers(i+aux1-1,j+aux2-1,k+aux3-1)[0]
                                  - centers(i,j,k)[0]
                                );

                            Aaux[index][1] =
                                weight *
                                (
                                    centers(i+aux1-1,j+aux2-1,k+aux3-1)[1]
                                  - centers(i,j,k)[1]
                                );

                            Aaux[index][2] =
                                weight *
                                (
                                    centers(i+aux1-1,j+aux2-1,k+aux3-1)[2]
                                  - centers(i,j,k)[2]
                                );
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
                        A[3*aux1 + aux2] +=
                            Aaux[aux3][aux1] * Aaux[aux3][aux2];
                    }

                    b[aux1] += Aaux[aux3][aux1] * baux[aux3];
                }
            }

            n(i,j,k) = A.inv() & b;
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
