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
        fvMsh_.msh().faceBoundaryType();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                bool aux = true;

                if ((fvMsh_[0].l() > 1) && (i != 1))
                {
                    aux =
                        i == 0
                      ? aux && faceType.left()
                      > domainBoundary::typeNumber
                      : aux && faceType.right()
                      > domainBoundary::typeNumber;
                }

                if ((fvMsh_[0].m() > 1) && (j != 1))
                {
                    aux =
                        j == 0
                      ? aux && faceType.bottom()
                      > domainBoundary::typeNumber
                      : aux && faceType.top()
                      > domainBoundary::typeNumber;
                }

                if ((fvMsh_[0].n() > 1) && (k != 1))
                {
                    aux =
                        k == 0
                      ? aux && faceType.aft()
                      > domainBoundary::typeNumber
                      : aux && faceType.fore()
                      > domainBoundary::typeNumber;
                }

                boundaryType_[i][j][k] = aux;
            }
        }
    }
}

LSGIR::LSGIR
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dict, alpha)
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

void LSGIR::correct()
{
    colocatedVectorField& n = *this;

    // Reconstruct the interface normal using the version of LSGIR presented in
    // Lopez (2022)

    const meshField<vector,colocated>& cc =
        fvMsh_.template metrics<colocated>().cellCenters();

    int index;
    scalar weight;

    forAllCells(n, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) > vof::threshold)
         && (alpha_(i,j,k) < 1 - vof::threshold)
        )
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
                                        cc(i,j,k)
                                      - cc(i+aux1-1,j+aux2-1,k+aux3-1)
                                    ),
                                    1.5
                                );

                            baux[index] =
                                weight *
                                (
                                    alpha_(i+aux1-1,j+aux2-1,k+aux3-1)
                                  - alpha_(i,j,k)
                                );

                            Aaux[index][0] =
                                weight *
                                (
                                    cc(i+aux1-1,j+aux2-1,k+aux3-1)[0]
                                  - cc(i,j,k)[0]
                                );

                            Aaux[index][1] =
                                weight *
                                (
                                    cc(i+aux1-1,j+aux2-1,k+aux3-1)[1]
                                  - cc(i,j,k)[1]
                                );

                            Aaux[index][2] =
                                weight *
                                (
                                    cc(i+aux1-1,j+aux2-1,k+aux3-1)[2]
                                  - cc(i,j,k)[2]
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

    n[0].correctBoundaryConditions();
}

}

}

}
