#include "SHF.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SHF, 0);
addToRunTimeSelectionTable(curvatureScheme, SHF, dictionary);

SHF::SHF
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    curvatureScheme(fvMsh, dict, normal, alpha)
{}

SHF::SHF(const SHF& s)
:
    curvatureScheme(s)
{}

SHF::~SHF()
{}

void SHF::correct()
{
    colocatedScalarField& kappa = *this;

    colocatedScalarField kc2("curvature2", fvMsh_);
    colocatedScalarField marker("marker", fvMsh_);

    const rectilinearMesh& rMsh = fvMsh_.msh().cast<rectilinearMesh>();

    const PartialList<scalar>& xSize = rMsh.localCellSizes()[0];
    const PartialList<scalar>& ySize = rMsh.localCellSizes()[1];
    const PartialList<scalar>& zSize = rMsh.localCellSizes()[2];

    const tensor T = rMsh.base();
    const tensor TtT = T & T.T();
    const scalar angleTol = Foam::cos(0.8);

    forAllCells(kappa, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) > vof::threshold)
         && (alpha_(i,j,k) < 1 - vof::threshold)
        )
        {

            int d1 = 0;
            int d2 = 1;
            int d3 = 2;
            scalar value;
            value = Foam::mag(normal_(i,j,k)[0]);
            int index1 = j;
            int index2 = k;

            scalar minH = 0;
            scalar maxH = 0;

            if (Foam::mag(normal_(i,j,k)[1]) > value)
            {
                d1 = 1;
                d2 = 0;
                value = Foam::mag(normal_(i,j,k)[1]);
                index1 = i;
            }

            if (Foam::mag(normal_(i,j,k)[2]) > value)
            {
                d1 = 2;
                d2 = 0;
                d3 = 1;
                value = Foam::mag(normal_(i,j,k)[2]);
                index1 = i;
                index2 = j;
            }

            scalar H[3][3] = {0};
            label tup = 0;
            label tdown = 0;
            label nsign = normal_(i,j,k)[d1] > 0 ? 1 : -1;

            scalar zeroSum = 0;
            scalar prevSum = 0;
            scalar newSum = 0;

            const scalar gamma = value > angleTol ? 0.0 : 0.2;

            if (d1 == 0)
            {
                for (int aux1 = -1; aux1 <= 1; aux1++)
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                        zeroSum += alpha_(i,j + aux1,k + aux2);

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i + aux3,j + aux1,k + aux2);

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                     && newSum != 0
                     && newSum != 9
                     && (i + aux3) <= kappa.I().right()
                    )
                    {
                        tup++;
                        prevSum = newSum;
                        maxH += xSize[i + aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i - aux3,j + aux1,k + aux2);

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                     && newSum != 0
                     && newSum != 9
                     && (i - aux3) >= (kappa.I().left() - 1)
                    )
                    {
                        tdown++;
                        prevSum = newSum;
                        minH += xSize[i - aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                if (nsign > 0)
                    minH = maxH;

                maxH  = minH + xSize[i];

                for (int aux1 = -1; aux1 <= 1; aux1++)
                {
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                    {
                        H[aux1 + 1][aux2 + 1] =
                            xSize[i] * alpha_(i, j + aux1, k + aux2);

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i+aux3,j+aux1,k+aux2)
                                  - alpha_(i+aux3-1,j+aux1,k+aux2)
                                ) > 0
                              ? xSize[i+aux3]*alpha_(i+aux3,j+aux1,k+aux2)
                              : nsign > 0 ? xSize[i+aux3] : 0;
                        }

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i-aux3,j+aux1,k+aux2)
                                  - alpha_(i-aux3+1,j+aux1,k+aux2)
                                ) < 0
                              ? xSize[i-aux3]*alpha_(i-aux3,j+aux1,k+aux2)
                              : nsign < 0 ? xSize[i-aux3] : 0;
                        }
                    }
                }
            }
            else if (d1 == 1)
            {
                for (int aux1 = -1; aux1 <= 1; aux1++)
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                        zeroSum += alpha_(i + aux1, j, k + aux2);

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i + aux1,j + aux3,k + aux2);

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                      && newSum != 0
                      && newSum != 9
                      && (j + aux3) <= kappa.I().top()
                    )
                    {
                        tup++;
                        prevSum = newSum;
                        maxH += ySize[j+aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i + aux1,j - aux3,k + aux2);

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                      && newSum != 0
                      && newSum != 9
                      && (j - aux3) >= (kappa.I().bottom() - 1)
                    )
                    {
                        tdown++;
                        prevSum = newSum;
                        minH += ySize[j-aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                if (nsign > 0)
                    minH = maxH;

                maxH  = minH + ySize[j];

                for (int aux1 = -1; aux1 <= 1; aux1++)
                {
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                    {
                        H[aux1 + 1][aux2 + 1] =
                            ySize[j] * alpha_(i + aux1, j, k + aux2);

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i+aux1,j+aux3,k+aux2)
                                  - alpha_(i+aux1,j+aux3-1,k+aux2)
                                ) > 0
                              ? ySize[j+aux3]*alpha_(i+aux1,j+aux3,k+aux2)
                              : nsign > 0 ? ySize[j+aux3] : 0;
                        }

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i+aux1,j-aux3,k+aux2)
                                  - alpha_(i+aux1,j-aux3+1,k+aux2)
                                ) < 0
                              ? ySize[j-aux3]*alpha_(i+aux1,j-aux3,k+aux2)
                              : nsign < 0 ? ySize[j-aux3] : 0;
                        }
                    }
                }
            }
            else
            {
                for (int aux1 = -1; aux1 <= 1; aux1++)
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                        zeroSum += alpha_(i + aux1,j + aux2,k);

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i + aux1,j + aux2,k + aux3);

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                      && newSum != 0
                      && newSum != 9
                      && (k + aux3) <= kappa.I().fore()
                    )
                    {
                        tup++;
                        prevSum = newSum;
                        maxH += zSize[k+aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                prevSum = zeroSum;

                for (int aux3 = 1; aux3 <=3; aux3++)
                {
                    newSum = 0;

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            newSum += alpha_(i + aux1,j + aux2,k - aux3);

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                     && newSum != 0
                     && newSum != 9
                     && (k - aux3) >= (kappa.I().aft() - 1)
                    )
                    {
                        tdown++;
                        prevSum = newSum;
                        minH += zSize[k-aux3];
                    }
                    else
                    {
                        break;
                    }
                }

                if (nsign > 0)
                    minH = maxH;

                maxH  = minH + zSize[k];

                for (int aux1 = -1; aux1 <= 1; aux1++)
                {
                    for (int aux2 = -1; aux2 <= 1; aux2++)
                    {
                        H[aux1 + 1][aux2 + 1] =
                            zSize[k] * alpha_(i + aux1, j + aux2, k);

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i+aux1,j+aux2,k+aux3)
                                  - alpha_(i+aux1,j+aux2,k+aux3-1)
                                ) > 0
                              ? zSize[k+aux3]*alpha_(i+aux1,j+aux2,k+aux3)
                              : nsign > 0 ? zSize[k+aux3] : 0;
                        }

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    alpha_(i+aux1,j+aux2,k-aux3)
                                  - alpha_(i+aux1,j+aux2,k-aux3+1)
                                ) < 0
                              ? zSize[k-aux3]*alpha_(i+aux1,j+aux2,k-aux3)
                              : nsign < 0 ? zSize[k-aux3] : 0;
                        }
                    }
                }
            }

            const PartialList<scalar>& xSizes = rMsh.localCellSizes()[d2];
            const PartialList<scalar>& ySizes = rMsh.localCellSizes()[d3];

            scalarList Hxx(3);
            scalarList Hyy(3);
            scalarList Hx(3);
            scalarList Hy(3);

            for (int aux1 = 0; aux1 < 3; aux1++)
            {
                Hxx[aux1] =
                    (1/xSizes[index1]) *
                    (
                        (H[2][aux1] - H[1][aux1])
                      / (0.5 * (xSizes[index1+1] + xSizes[index1]))
                      - (H[1][aux1] - H[0][aux1])
                      / (0.5 * (xSizes[index1] + xSizes[index1-1]))
                    );

                Hyy[aux1] =
                    (1/ySizes[index2]) *
                    (
                        (H[aux1][2] - H[aux1][1])
                      / (0.5 * (ySizes[index2+1] + ySizes[index2]))
                      - (H[aux1][1] - H[aux1][0])
                      / (0.5 * (ySizes[index2] + ySizes[index2-1]))
                    );

                Hx[aux1] =
                    H[0][aux1] *
                    (
                        (-xSizes[index1]-xSizes[index1+1])
                      / (
                            (xSizes[index1]+xSizes[index1-1])
                          * (
                                0.5*(xSizes[index1-1]+xSizes[index1+1])
                              + xSizes[index1]
                            )
                        )
                    )
                  + H[1][aux1] *
                    (
                        2*(xSizes[index1+1]-xSizes[index1-1])
                      / (
                            (xSizes[index1]+xSizes[index1+1])
                          * (xSizes[index1]+xSizes[index1-1])
                        )
                    )
                  + H[2][aux1] *
                    (
                        (xSizes[index1]+xSizes[index1-1])
                      / (
                            (xSizes[index1]+xSizes[index1+1])
                          * (
                                0.5*(xSizes[index1+1]+xSizes[index1-1])
                              + xSizes[index1]
                            )
                        )
                    );

                Hy[aux1] =
                    H[aux1][0] *
                    (
                        (-ySizes[index2]-ySizes[index2+1])
                      / (
                            (ySizes[index2]+ySizes[index2-1])
                          * (
                                0.5*(ySizes[index2-1]+ySizes[index2+1])
                              + ySizes[index2]
                            )
                        )
                    )
                  + H[aux1][1] *
                    (
                        2*(ySizes[index2+1]-ySizes[index2-1])
                      / (
                            (ySizes[index2]+ySizes[index2+1])
                          * (ySizes[index2]+ySizes[index2-1])
                        )
                    )
                    + H[aux1][2] *
                    (
                        (ySizes[index2]+ySizes[index2-1])
                      / (
                            (ySizes[index2]+ySizes[index2+1])
                          * (
                                0.5*(ySizes[index2+1]+ySizes[index2-1])
                              + ySizes[index2]
                            )
                        )
                    );
            }

            scalar Hxy =
                Hx[0] *
                (
                    (-ySizes[index2]-ySizes[index2+1])
                  / (
                        (ySizes[index2]+ySizes[index2-1])
                      * (0.5*(ySizes[index2-1]+ySizes[index2+1])+ySizes[index2])
                    )
                )
              + Hx[1] *
                (
                    2*(ySizes[index2+1]-ySizes[index2-1])
                  / (
                        (ySizes[index2]+ySizes[index2+1])
                      * (ySizes[index2]+ySizes[index2-1])
                    )
                )
              + Hx[2] *
                (
                    (ySizes[index2]+ySizes[index2-1])
                  / (
                        (ySizes[index2]+ySizes[index2+1])
                      * (0.5*(ySizes[index2+1]+ySizes[index2-1])+ySizes[index2])
                    )
                );

            if ((minH <= H[1][1]) && (maxH >= H[1][1]))
            {
                scalar dx = (gamma * (Hx[0] + Hx[2]) + Hx[1]) / (1 + 2 * gamma);
                scalar dy = (gamma * (Hy[0] + Hy[2]) + Hy[1]) / (1 + 2 * gamma);

                scalar dxx =
                    (gamma * (Hxx[0] + Hxx[2]) + Hxx[1]) / (1 + 2 * gamma);
                scalar dyy =
                    (gamma * (Hyy[0] + Hyy[2]) + Hyy[1]) / (1 + 2 * gamma);

                vector normal;
                normal[d2] = -dx;
                normal[d3] = -dy;
                normal[d1] = 1;
                normal = T.inv() & normal;
                normal /= Foam::mag(normal);

                vector dfdx = Zero;
                vector dfdy = Zero;
                vector dfdxx = Zero;
                vector dfdxy = Zero;
                vector dfdyy = Zero;

                dfdx[d1] = dx;
                dfdy[d1] = dy;
                dfdx[d2] = 1;
                dfdy[d3] = 1;
                dfdxx[d1] = dxx;
                dfdxy[d1] = Hxy;
                dfdyy[d1] = dyy;

                scalar kE = dfdx & (TtT & dfdx);
                scalar kF = dfdx & (TtT & dfdy);
                scalar kG = dfdy & (TtT & dfdy);

                scalar kL = dfdxx & (T & normal);
                scalar kM = dfdxy & (T & normal);
                scalar kN = dfdyy & (T & normal);


                kc2(i,j,k) =
                    (
                        kG * kL - 2 * kF * kM + kE * kN
                    )
                  / (
                        kE * kG - Foam::sqr(kF)
                    );

                marker(i,j,k) = 1;
            }
            else
            {
                kc2(i,j,k) = 0;
                marker(i,j,k) = 0;
            }
        }
        else
        {
            kappa(i,j,k) = 0;
            marker(i,j,k) = 0;
        }
    }

    forAllCells(kappa, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) > vof::threshold)
         && (alpha_(i,j,k) < 1 - vof::threshold)
        )
        {
            int count = 0;
            kappa(i,j,k) = 0;

            for (int aux1 = -1; aux1 <= 1; aux1++)
            {
                for (int aux2 = -1; aux2 <= 1; aux2++)
                {
                    if
                    (
                        Foam::mag(marker(i+aux1,j+aux2,k) - 1) < 1e-12
                    )
                    {
                        count++;
                        kappa(i,j,k) += kc2(i+aux1,j+aux2,k);
                    }
                }
            }

            kappa(i,j,k) /= scalar(count);
        }
    }
}

}

}

}
