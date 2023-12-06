#include "SHF.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"
#include "cellDataExchange.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SHF, 0);
addToRunTimeSelectionTable(curvatureScheme, SHF, dictionary);

void SHF::setCommunicationIndices()
{
    const faceLabel& faceType = fvMsh_.msh().faceBoundaryType();
    const faceLabel I = fvMsh_.template I<colocated>(0, 0);

    communicationIndices_.clear();
    communicationIndices_.setSize(6);

    if (faceType.left() > 1)
        for (int i = I.left() - 3; i < I.left() - 1; i++)
            for (int j = I.bottom() - 1; j < I.top() + 1; j++)
                for (int k = I.aft() - 1; k < I.fore() + 1; k++)
                    communicationIndices_[0].append(labelVector(i,j,k));

    if (faceType.right() > 1)
        for (int i = I.right() + 1; i < I.right() + 3; i++)
            for (int j = I.bottom() - 1; j < I.top() + 1; j++)
                for (int k = I.aft() - 1; k < I.fore() + 1; k++)
                    communicationIndices_[1].append(labelVector(i,j,k));

    if (faceType.bottom() > 1)
        for (int i = I.left() - 1; i < I.right(); i++)
            for (int j = I.bottom() - 3; j < I.bottom() - 1; j++)
                for (int k = I.aft() - 1; k < I.fore() + 1; k++)
                    communicationIndices_[2].append(labelVector(i,j,k));

    if (faceType.top() > 1)
        for (int i = I.left() - 1; i < I.right(); i++)
            for (int j = I.top() + 1; j < I.top() + 3; j++)
                for (int k = I.aft() - 1; k < I.fore() + 1; k++)
                    communicationIndices_[3].append(labelVector(i,j,k));

    if (faceType.aft() > 1)
        for (int i = I.left() - 1; i < I.right(); i++)
            for (int j = I.bottom() - 1; j < I.top(); j++)
                for (int k = I.aft() - 3; k < I.aft() - 1; k++)
                    communicationIndices_[4].append(labelVector(i,j,k));

    if (faceType.fore() > 1)
        for (int i = I.left() - 1; i < I.right(); i++)
            for (int j = I.bottom() - 1; j < I.top(); j++)
                for (int k = I.fore() + 1; k < I.fore() + 3; k++)
                    communicationIndices_[5].append(labelVector(i,j,k));

}

void SHF::setMaxKappa()
{
    const rectilinearMesh& rMsh = fvMsh_.msh().cast<rectilinearMesh>();

    const PartialList<scalar>& xSize = rMsh.localCellSizes()[0];
    const PartialList<scalar>& ySize = rMsh.localCellSizes()[1];
    const PartialList<scalar>& zSize = rMsh.localCellSizes()[2];

    maxKappa_ = Foam::max(Foam::max(1.0/xSize[0],1.0/ySize[0]),1.0/zSize[0]);
    reduce(maxKappa_, minOp<scalar>());
}

SHF::SHF
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    curvatureScheme(fvMsh, dict, normal, alpha)
{
    setCommunicationIndices();
    setMaxKappa();
}

SHF::SHF(const SHF& s)
:
    curvatureScheme(s),
    communicationIndices_(s.communicationIndices_),
    maxKappa_(s.maxKappa_)
{}

SHF::~SHF()
{}

void SHF::correct()
{
    colocatedScalarField& kappa = *this;

    const faceLabel& faceType = fvMsh_.msh().faceBoundaryType();
    const labelVector& N = fvMsh_.template N<colocated>(0, 0);

    colocatedScalarField marker("marker", fvMsh_);

    const rectilinearMesh& rMsh = fvMsh_.msh().cast<rectilinearMesh>();

    const PartialList<scalar>& xSize = rMsh.localCellSizes()[0];
    const PartialList<scalar>& ySize = rMsh.localCellSizes()[1];
    const PartialList<scalar>& zSize = rMsh.localCellSizes()[2];

    const scalar angleTol = Foam::cos(0.8);

    commKappa_.clear();

    for (int i = 0; i < 6; i++)
    {
        cellDataExchange<colocated> exchange(communicationIndices_[i], fvMsh_, 0, 0);
        List<scalar> data(move(exchange(kappa)));
        commKappa_.append(data);
    }

    forAllCells(kappa, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) >     vof::threshold)
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

                    if
                    (
                        (i + aux3) > kappa.I().right()
                     && faceType.right() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((i + aux3) <= kappa.I().right())
                            {
                                newSum += alpha_(i + aux3,j + aux1,k + aux2);
                            }
                            else
                            {
                                int index = (i + aux3 - kappa.I().right() - 1)
                                          * (N[1] + 2) * (N[2] + 2)
                                          + (j + aux1) * (N[2] + 2) + k + aux2;
                                newSum += commKappa_[1][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                     && newSum != 0
                     && newSum != 9
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

                    if
                    (
                        (i - aux3) < (kappa.I().left() - 1)
                     && faceType.left() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((i - aux3) >= (kappa.I().left() - 1))
                            {
                                newSum += alpha_(i - aux3,j + aux1,k + aux2);
                            }
                            else
                            {
                                int index = (i - aux3 - kappa.I().left() + 3)
                                          * (N[1] + 2) * (N[2] + 2)
                                          + (j + aux1) * (N[2] + 2) + k + aux2;
                                newSum += commKappa_[0][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                     && newSum != 0
                     && newSum != 9
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
                        scalar value = 0;
                        scalar valuePrev = alpha_(i, j + aux1, k + aux2);
                        H[aux1 + 1][aux2 + 1] =
                            xSize[i] * valuePrev;

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            if ((i + aux3) <= kappa.I().right())
                            {
                                value = alpha_(i + aux3,j + aux1,k + aux2);
                            }
                            else
                            {
                                int index = (i + aux3 - kappa.I().right() - 1)
                                          * (N[1] + 2) * (N[2] + 2)
                                          + (j + aux1) * (N[2] + 2) + k + aux2;
                                value = commKappa_[1][index];
                            }

                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) > 0
                              ? xSize[i+aux3]*value
                              : nsign > 0 ? xSize[i+aux3] : 0;

                            valuePrev = value;
                        }

                        valuePrev = alpha_(i, j + aux1, k + aux2);

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            if ((i - aux3) >= (kappa.I().left() - 1))
                            {
                                value = alpha_(i - aux3,j + aux1,k + aux2);
                            }
                            else
                            {
                                int index = (i - aux3 - kappa.I().left() + 3)
                                          * (N[1] + 2) * (N[2] + 2)
                                          + (j + aux1) * (N[2] + 2) + k + aux2;
                                value = commKappa_[0][index];
                            }

                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) < 0
                              ? xSize[i-aux3]*value
                              : nsign < 0 ? xSize[i-aux3] : 0;

                            valuePrev = value;
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

                    if
                    (
                        (j + aux3) > kappa.I().top()
                     && faceType.top() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((j + aux3) <= kappa.I().top())
                            {
                                newSum += alpha_(i + aux1,j + aux3,k + aux2);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[2] + 2)
                                          + (j + aux3 - kappa.I().top() - 1)
                                          * (N[2] + 2) + k + aux2;
                                newSum += commKappa_[3][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                      && newSum != 0
                      && newSum != 9
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

                    if
                    (
                        (j - aux3) < (kappa.I().bottom() - 1)
                     && faceType.bottom() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((j - aux3) >= (kappa.I().bottom() - 1))
                            {
                                newSum += alpha_(i + aux1,j - aux3,k + aux2);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[2] + 2)
                                          + (j - aux3 - kappa.I().bottom() + 3)
                                          * (N[2] + 2) + k + aux2;
                                newSum += commKappa_[2][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                      && newSum != 0
                      && newSum != 9
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
                        scalar value = 0;
                        scalar valuePrev = alpha_(i + aux1, j, k + aux2);

                        H[aux1 + 1][aux2 + 1] =
                            ySize[j] * valuePrev;

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            if ((j + aux3) <= kappa.I().top())
                            {
                                value = alpha_(i + aux1,j + aux3,k + aux2);
                            }
                            else
                            {

                                int index = (i + aux1) * 2 * (N[2] + 2)
                                          + (j + aux3 - kappa.I().top() - 1)
                                          * (N[2] + 2) + k + aux2;
                                value = commKappa_[3][index];
                            }

                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) > 0
                              ? ySize[j+aux3]*value
                              : nsign > 0 ? ySize[j+aux3] : 0;

                            valuePrev = value;
                        }

                        valuePrev = alpha_(i + aux1, j, k + aux2);

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            if ((j - aux3) >= (kappa.I().bottom() - 1))
                            {
                                value = alpha_(i + aux1,j - aux3,k + aux2);
                            }
                            else
                            {

                                int index = (i + aux1) * 2 * (N[2] + 2)
                                          + (j - aux3 - kappa.I().bottom() + 3)
                                          * (N[2] + 2) + k + aux2;
                                value = commKappa_[2][index];
                            }
                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) < 0
                              ? ySize[j-aux3]*value
                              : nsign < 0 ? ySize[j-aux3] : 0;

                            valuePrev = value;
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

                    if
                    (
                        (k + aux3) > kappa.I().fore()
                     && faceType.fore() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((k + aux3) <= kappa.I().fore())
                            {
                                newSum += alpha_(i + aux1,j + aux2,k + aux3);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[1] + 2)
                                          + (j + aux2) * 2
                                          + k + aux3 - kappa.I().fore() - 1;
                                newSum += commKappa_[5][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) >= 0
                      && newSum != 0
                      && newSum != 9
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

                    if
                    (
                        (k - aux3) < (kappa.I().aft() - 1)
                     && faceType.fore() < 2
                    )
                    {
                        break;
                    }

                    for (int aux1 = -1; aux1 <= 1; aux1++)
                        for (int aux2 = -1; aux2 <= 1; aux2++)
                            if ((k - aux3) >= (kappa.I().aft() - 1))
                            {
                                newSum += alpha_(i + aux1,j + aux2,k - aux3);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[1] + 2)
                                          + (j + aux2) * 2
                                          + k - aux3 - kappa.I().aft() + 3;
                                newSum += commKappa_[4][index];
                            }

                    if
                    (
                        nsign * (newSum - prevSum) <= 0
                     && newSum != 0
                     && newSum != 9
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
                        scalar value = 0;
                        scalar valuePrev = alpha_(i + aux1, j + aux2, k);

                        H[aux1 + 1][aux2 + 1] =
                            zSize[k] * valuePrev;

                        for (int aux3 = 1; aux3 <= tup; aux3++)
                        {
                            if ((k + aux3) <= kappa.I().fore())
                            {
                                value = alpha_(i + aux1,j + aux2,k + aux3);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[1] + 2)
                                          + (j + aux2) * 2
                                          + k + aux3 - kappa.I().fore() - 1;
                                value = commKappa_[5][index];
                            }

                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) > 0
                              ? zSize[k+aux3]*value
                              : nsign > 0 ? zSize[k+aux3] : 0;

                            valuePrev = value;
                        }

                        valuePrev = alpha_(i + aux1, j + aux2, k);

                        for (int aux3 = 1; aux3 <= tdown; aux3++)
                        {
                            if ((k - aux3) >= (kappa.I().aft() - 1))
                            {
                                value = alpha_(i+aux1,j+aux2,k-aux3);
                            }
                            else
                            {
                                int index = (i + aux1) * 2 * (N[1] + 2)
                                          + (j + aux2) * 2
                                          + k - aux3 - kappa.I().aft() + 3;
                                value = commKappa_[4][index];
                            }

                            H[aux1 + 1][aux2 + 1] +=
                                nsign
                              * (
                                    value
                                  - valuePrev
                                ) < 0
                              ? zSize[k-aux3]*value
                              : nsign < 0 ? zSize[k-aux3] : 0;

                            valuePrev = value;
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

                kappa(i,j,k) = -
                    (
                        dxx * (1 + Foam::sqr(dy))
                      + dyy * (1 + Foam::sqr(dx))
                      - 2 * Hxy * dx * dy
                    )
                  / (
                        Foam::pow(1 + Foam::sqr(dx) + Foam::sqr(dy), 1.5)
                    );

                if (Foam::mag(kappa(i,j,k)) > maxKappa_)
                {
                    kappa(i,j,k) = maxKappa_ * Foam::sign(kappa(i,j,k));
                }

                marker(i,j,k) = 1;
            }
            else
            {
                marker(i,j,k) = 0;
            }
        }
        else
        {
            kappa(i,j,k) = 0;
            marker(i,j,k) = 0;
        }
    }

    marker.correctCommBoundaryConditions();
    kappa.correctCommBoundaryConditions();

    forAllCells(kappa, i, j, k)
    {
        if
        (
            (alpha_(i,j,k) >     vof::threshold)
         && (alpha_(i,j,k) < 1 - vof::threshold)
         && (Foam::mag(marker(i,j,k)) < 1e-12)
        )
        {
            int count = 0;
            kappa(i,j,k) = 0;

            for (int aux1 = -1; aux1 <= 1; aux1++)
            {
                for (int aux2 = -1; aux2 <= 1; aux2++)
                {
                    for (int aux3 = -1; aux3 <= 1; aux3++)
                    {
                        if
                        (
                            Foam::mag(marker(i+aux1,j+aux2,k+aux3) - 1) < 1e-12
                        )
                        {
                            count++;
                            kappa(i,j,k) += kappa(i+aux1,j+aux2,k+aux3);
                        }
                    }
                }
            }

            if (count > 0)
            {
                kappa(i,j,k) /= scalar(count);
            }
        }
    }

    kappa.correctBoundaryConditions();
}

}

}

}
