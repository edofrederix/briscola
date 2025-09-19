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

void SHF::setCommunication()
{
    const faceLabel& boundaryType = fvMsh_.msh().faceBoundaryType();
    const faceLabel I = fvMsh_.template I<colocated>();

    List<labelVector> indices(0);

    exchangeShapes_.setSize(6, Zero);

    label cursor = 0;

    for (int f = 0; f < 6; f++)
    {
        const label d = f/2;

        // On parallel/periodic boundaries only

        if (boundaryType[f] > patchBoundary::typeNumber)
        {
            faceLabel J(I.lower() - unitXYZ, I.upper() + unitXYZ);

            if (f%2 == 0)
            {
                J[d*2  ] = I[f] - 3;
                J[d*2+1] = I[f] - 1;
            }
            else
            {
                J[d*2  ] = I[f] + 1;
                J[d*2+1] = I[f] + 3;
            }

            exchangeShapes_[f] = J.upper() - J.lower();

            indices.setSize
            (
                indices.size()
              + cmptProduct(exchangeShapes_[f])
            );

            for (int i = J.left(); i < J.right(); i++)
                for (int j = J.bottom(); j < J.top(); j++)
                    for (int k = J.aft(); k < J.fore(); k++)
                        indices[cursor++] = labelVector(i,j,k);
        }
    }

    exchangePtr_.reset
    (
        new cellDataExchange<colocated>(indices, fvMsh_)
    );
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
    setCommunication();
}

SHF::SHF(const SHF& s)
:
    curvatureScheme(s)
{
    setCommunication();
}

SHF::~SHF()
{}

void SHF::correct()
{
    colocatedScalarField& kappa = *this;
    kappa = Zero;

    const faceLabel& boundaryType = fvMsh_.msh().faceBoundaryType();
    const labelVector& N = fvMsh_.template N<colocated>();
    const faceLabel& I = fvMsh_.template I<colocated>();

    colocatedLabelField marker("marker", fvMsh_);
    marker = Zero;

    const rectilinearMesh& rMsh = fvMsh_.msh().cast<rectilinearMesh>();

    const PartialList<scalar>& xSize = rMsh.localCellSizes()[0];
    const PartialList<scalar>& ySize = rMsh.localCellSizes()[1];
    const PartialList<scalar>& zSize = rMsh.localCellSizes()[2];

    const labelVector& globalStart = rMsh.globalStarts()[Pstream::myProcNo()];

    // Make a padded block for alpha. The default block contains the ghost
    // cells. This block is then padded to contain boundary exchange cells.

    labelVector M(N + 2*unitXYZ);
    labelVector off(unitXYZ);

    for (int f = 0; f < 6; f++)
        if (exchangeShapes_[f][f/2] > 1)
            M[f/2] += exchangeShapes_[f][f/2];

    for (int d = 0; d < 3; d++)
        if (exchangeShapes_[d*2][d] > 1)
            off[d] += exchangeShapes_[d*2][d];

    scalarBlock alpha(M, Zero);

    // Copy data including ghost cells

    for (int i = -1; i < N.x() + 1; i++)
        for (int j = -1; j < N.y() + 1; j++)
            for (int k = -1; k < N.z() + 1; k++)
                alpha(labelVector(i,j,k) + off) = alpha_(i,j,k);

    // Fill with neighbor exchange data

    if (exchangePtr_.valid())
    {
        scalarList data(move(exchangePtr_()(alpha_)));

        for (int i = 0; i < data.size(); i++)
            alpha(exchangePtr_->indices()[i] + off) = data[i];
    }

    // Compute the curvature

    forAllCells(alpha_, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if
        (
            (alpha_(ijk) >     vof::threshold)
         && (alpha_(ijk) < 1 - vof::threshold)
        )
        {
            const vector m(cmptMag(normal_(ijk)));

            // Determine the primary interfacial orientation

            const label p =
                m.x() == cmptMax(m) ? 0
              : m.y() == cmptMax(m) ? 1
              :                       2;

            // Rotated indices based on primary direction

            labelVector d;

            for (int ii = 0; ii < 3; ii++)
                d[ii] = (ii + p) % 3;

            // Rotated base

            const labelTensor D(units[d[0]], units[d[1]], units[d[2]]);

            // Rotated face labels

            faceLabel f;

            for (int ii = 0; ii < 3; ii++)
            {
                f[ii*2    ] = d[ii]*2;
                f[ii*2 + 1] = d[ii]*2 + 1;
            }

            const PartialList<scalar>& gSizes = rMsh.globalCellSizes()[p];

            scalar minH = 0.0;
            scalar maxH = 0.0;

            List<scalarList> H(3, scalarList(3, 0.0));

            label upper = 0.0;
            label lower = 0.0;

            label s = Foam::sign(normal_(ijk)[p]);

            const scalar gamma = m[p] > Foam::cos(0.8) ? 0.0 : 0.2;

            scalar sum = 0;

            for (int b = -1; b <= 1; b++)
                for (int c = -1; c <= 1; c++)
                    sum += alpha(ijk + b*D.y() + c*D.z() + off);

            scalar sumPrev = sum;
            scalar sumNew;

            // Rotated cell indices

            const label ii = ijk[d.x()];
            const label jj = ijk[d.y()];
            const label kk = ijk[d.z()];

            // Look in positive primary direction

            for (int a = 1; a <= 3; a++)
            {
                sumNew = 0.0;

                // Do not go beyond domain boundary ghost cells

                if
                (
                    ii + a > I[f.right()]
                 && boundaryType[f.right()] < parallelBoundary::typeNumber
                )
                {
                    break;
                }

                for (int b = -1; b <= 1; b++)
                    for (int c = -1; c <= 1; c++)
                        sumNew +=
                            alpha(ijk + a*D.x() + b*D.y() + c*D.z() + off);

                if
                (
                    s*(sumNew - sumPrev) >= 0.0
                 && sumNew >       1.0e-12
                 && sumNew < 9.0 - 1.0e-12
                )
                {
                    upper++;
                    sumPrev = sumNew;
                    maxH += gSizes[globalStart[p] + ii + a];
                }
                else
                {
                    break;
                }
            }

            sumPrev = sum;

            // Look in negative primary direction

            for (int a = 1; a <= 3; a++)
            {
                sumNew = 0.0;

                // Do not go beyond domain boundary ghost cells

                if
                (
                    ii - a < (I[f.left()] - 1)
                 && boundaryType[f.left()] < parallelBoundary::typeNumber
                )
                {
                    break;
                }

                for (int b = -1; b <= 1; b++)
                    for (int c = -1; c <= 1; c++)
                        sumNew +=
                            alpha(ijk - a*D.x() + b*D.y() + c*D.z() + off);

                if
                (
                    s*(sumNew - sumPrev) <= 0.0
                 && sumNew >       1.0e-12
                 && sumNew < 9.0 - 1.0e-12
                )
                {
                    lower++;
                    sumPrev = sumNew;
                    minH += gSizes[globalStart[p] + ii - a];
                }
                else
                {
                    break;
                }
            }

            if (s > 0)
                minH = maxH;

            maxH  = minH + gSizes[globalStart[p] + ii];

            for (int b = -1; b <= 1; b++)
            {
                for (int c = -1; c <= 1; c++)
                {
                    scalar value = 0;
                    scalar valuePrev = alpha(ijk + b*D.y() + c*D.z() + off);

                    H[b+1][c+1] = gSizes[globalStart[p] + ii]*valuePrev;

                    // Upper

                    for (int a = 1; a <= upper; a++)
                    {
                        value = alpha(ijk + a*D.x() + b*D.y() + c*D.z() + off);

                        H[b+1][c+1] +=
                            s*(value - valuePrev) > 0.0
                          ? gSizes[globalStart[p] + ii + a]*value
                          : s > 0.0
                          ? gSizes[globalStart[p] + ii + a]
                          : 0.0;

                        valuePrev = value;
                    }

                    valuePrev = alpha(ijk + b*D.y() + c*D.z() + off);

                    // Lower

                    for (int a = 1; a <= lower; a++)
                    {
                        value = alpha(ijk - a*D.x() + b*D.y() + c*D.z() + off);

                        H[b+1][c+1] +=
                            s*(value - valuePrev) < 0.0
                          ? gSizes[globalStart[p] + ii - a]*value
                          : s < 0.0
                          ? gSizes[globalStart[p] + ii - a]
                          : 0.0;

                        valuePrev = value;
                    }
                }
            }

            const PartialList<scalar>& dy = rMsh.localCellSizes()[d.y()];
            const PartialList<scalar>& dz = rMsh.localCellSizes()[d.z()];

            scalarList Hyy(3);
            scalarList Hzz(3);
            scalarList Hy(3);
            scalarList Hz(3);

            for (int b = 0; b < 3; b++)
            {
                Hyy[b] =
                    1.0/dy[jj]
                  * (
                        (H[2][b] - H[1][b])/(0.5*(dy[jj+1] + dy[jj  ]))
                      - (H[1][b] - H[0][b])/(0.5*(dy[jj  ] + dy[jj-1]))
                    );

                Hzz[b] =
                    1.0/dz[kk]
                  * (
                        (H[b][2] - H[b][1])/(0.5*(dz[kk+1] + dz[kk  ]))
                      - (H[b][1] - H[b][0])/(0.5*(dz[kk  ] + dz[kk-1]))
                    );

                Hy[b] =
                  - H[0][b]*(dy[jj] + dy[jj+1])
                  / ((dy[jj] + dy[jj-1])*(0.5*(dy[jj-1] + dy[jj+1]) + dy[jj]))
                  + H[1][b]*2.0*(dy[jj+1] - dy[jj-1])
                  / ((dy[jj] + dy[jj+1])*(dy[jj] + dy[jj-1]))
                  + H[2][b]*(dy[jj] + dy[jj-1])
                  / ((dy[jj] + dy[jj+1])*(0.5*(dy[jj+1] + dy[jj-1]) + dy[jj]));

                Hz[b] =
                  - H[b][0]*(dz[kk] + dz[kk+1])
                  / ((dz[kk] + dz[kk-1])*(0.5*(dz[kk-1] + dz[kk+1]) + dz[kk]))
                  + H[b][1]*2.0*(dz[kk+1] - dz[kk-1])
                  / ((dz[kk] + dz[kk+1])*(dz[kk] + dz[kk-1]))
                  + H[b][2]*(dz[kk] + dz[kk-1])
                  / ((dz[kk] + dz[kk+1])*(0.5*(dz[kk+1] + dz[kk-1]) + dz[kk]));
            }

            const scalar Hyz =
              - Hy[0]*(dz[kk] + dz[kk+1])
              / ((dz[kk] + dz[kk-1])*(0.5*(dz[kk-1] + dz[kk+1]) + dz[kk]))

              + Hy[1]*2.0*(dz[kk+1] - dz[kk-1])
              / ((dz[kk] + dz[kk+1])*(dz[kk] + dz[kk-1]))

              + Hy[2]*(dz[kk] + dz[kk-1])
              / ((dz[kk] + dz[kk+1])*(0.5*(dz[kk+1] + dz[kk-1]) + dz[kk]));

            if (minH <= H[1][1] && maxH >= H[1][1])
            {
                const scalar Dy =
                    (gamma*(Hy[0] + Hy[2]) + Hy[1])/(1.0 + 2.0*gamma);
                const scalar Dz =
                    (gamma*(Hz[0] + Hz[2]) + Hz[1])/(1.0 + 2.0*gamma);

                const scalar Dyy =
                     (gamma*(Hyy[0] + Hyy[2]) + Hyy[1])/(1.0 + 2.0*gamma);
                const scalar Dzz =
                     (gamma*(Hzz[0] + Hzz[2]) + Hzz[1])/(1.0 + 2.0*gamma);

                kappa(ijk) =
                  - (
                        Dyy*(1.0 + Foam::sqr(Dz))
                      + Dzz*(1.0 + Foam::sqr(Dy))
                      - 2.0*Hyz*Dy*Dz
                    )
                  / Foam::pow(1.0 + Foam::sqr(Dy) + Foam::sqr(Dz), 1.5);

                // Limit kappa

                const scalar kappaMax =
                    1.0/cmptMin(vector(xSize[i], ySize[j], zSize[k]));

                if (Foam::mag(kappa(ijk)) > kappaMax)
                    kappa(ijk) = kappaMax * Foam::sign(kappa(ijk));

                // Set marker

                marker(ijk) = 1;
            }
        }
    }

    marker.correctCommsBoundaryConditions();
    kappa.correctCommsBoundaryConditions();

    // Set kappa in unmarked interfacial cells to be the average of its marked
    // neighbors

    forAllCells(alpha_, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if
        (
            (alpha_(ijk) >     vof::threshold)
         && (alpha_(ijk) < 1 - vof::threshold)
         && !marker(ijk)
        )
        {
            int count = 0;

            kappa(ijk) = 0.0;

            for (int b = -1; b <= 1; b++)
            {
                for (int c = -1; c <= 1; c++)
                {
                    for (int a = -1; a <= 1; a++)
                    {
                        const labelVector abc(a,b,c);

                        if (marker(ijk + abc))
                        {
                            kappa(ijk) += kappa(ijk + abc);
                            count++;
                        }
                    }
                }
            }

            if (count > 0)
                kappa(ijk) /= scalar(count);
        }
    }

    kappa.correctBoundaryConditions();
}

}

}

}
