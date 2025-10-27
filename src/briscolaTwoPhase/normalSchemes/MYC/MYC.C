#include "MYC.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"
#include "exSchemes.H"
#include "rectilinearMesh.H"

#include <cmath>

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(MYC, 0);
addToRunTimeSelectionTable(normalScheme, MYC, dictionary);

MYC::MYC
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    normalScheme(fvMsh, dict, alpha)
{}

MYC::MYC(const MYC& s)
:
    normalScheme(s)
{}

MYC::~MYC()
{}

void MYC::correct()
{
    colocatedVectorField& n = *this;

    // Reconstruct the normal interface using the Mixed Youngs method of Aulisa
    // et al. (2007).
    //
    // The method has been generalized for non uniform meshes using three points
    // Lagrange derivative approximations and for skewed meshes by transforming
    // the coordinate system into a perpendicular one.

    const rectilinearMesh& rMsh = fvMsh_.msh().cast<rectilinearMesh>();
    tensor T = rMsh.base();

    const PartialList<scalar>& size0 = rMsh.localCellSizes()[0];
    const PartialList<scalar>& size1 = rMsh.localCellSizes()[1];
    const PartialList<scalar>& size2 = rMsh.localCellSizes()[2];

    forAllCells(n, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if
        (
            (alpha_(ijk) > vof::threshold)
         && (alpha_(ijk) < 1 - vof::threshold)
        )
        {
            scalar m1, m2, m3;
            scalar t0, t1, t2;

            List<scalarList> m(4, scalarList(3));

            // Write the plane as:
            // sgn(mx)*X = my*Y + mz*Z + alpha*m00*X = m01*Y + m02*Z + alpha

            m1 =
                size0[i-1]
              * (
                    size1[j]  *size2[k]  *alpha_(i-1,j,  k  )
                  + size1[j+1]*size2[k]  *alpha_(i-1,j+1,k  )
                  + size1[j-1]*size2[k]  *alpha_(i-1,j-1,k  )
                  + size1[j]  *size2[k+1]*alpha_(i-1,j,  k+1)
                  + size1[j]  *size2[k-1]*alpha_(i-1,j,  k-1)
                );

            m2 =
                size0[i+1]
              * (
                    size1[j]  *size2[k]  *alpha_(i+1,j,  k  )
                  + size1[j+1]*size2[k]  *alpha_(i+1,j+1,k  )
                  + size1[j-1]*size2[k]  *alpha_(i+1,j-1,k  )
                  + size1[j]  *size2[k+1]*alpha_(i+1,j,  k+1)
                  + size1[j]  *size2[k-1]*alpha_(i+1,j,  k-1)
                );

            m[0][0] = m1 > m2 ? -1.0 : 1.0;

            m1 =
                size0[i-1]*alpha_(i-1,j-1,k)
              + size0[i+1]*alpha_(i+1,j-1,k)
              + size0[i  ]*alpha_(i,  j-1,k);

            m2 =
                size0[i-1]*alpha_(i-1,j,k)
              + size0[i+1]*alpha_(i+1,j,k)
              + size0[i  ]*alpha_(i,  j,k);

            m3 =
                size0[i-1]*alpha_(i-1,j+1,k)
              + size0[i+1]*alpha_(i+1,j+1,k)
              + size0[i  ]*alpha_(i,  j+1,k);

            m[0][1] =
                m1
              * (
                  - (size1[j] + size1[j+1])
                  / (
                        (size1[j] + size1[j-1])
                      * (0.5*(size1[j-1] + size1[j+1]) + size1[j])
                    )
                )
              + m2
              * (
                    2.0*(size1[j+1] - size1[j-1])
                  / (
                        (size1[j] + size1[j+1])
                      * (size1[j] + size1[j-1])
                    )
                )
              + m3
              * (
                    (size1[j] + size1[j-1])
                  / (
                        (size1[j] + size1[j+1])
                      * (0.5*(size1[j+1] + size1[j-1]) + size1[j])
                    )
                );

            m1 =
                size0[i-1]*alpha_(i-1,j,k-1)
              + size0[i+1]*alpha_(i+1,j,k-1)
              + size0[i  ]*alpha_(i,  j,k-1);

            m2 =
                size0[i-1]*alpha_(i-1,j,k)
              + size0[i+1]*alpha_(i+1,j,k)
              + size0[i  ]*alpha_(i,j,k);

            m3 =
                size0[i-1]*alpha_(i-1,j,k+1)
              + size0[i+1]*alpha_(i+1,j,k+1)
              + size0[i  ]*alpha_(i,  j,k+1);

            m[0][2] =
                m1
              * (
                  - (size2[k] + size2[k+1])
                  / (
                        (size2[k] + size2[k-1])
                      * (0.5*(size2[k-1] + size2[k+1]) + size2[k])
                    )
                )
              + m2
              * (
                    2.0*(size2[k+1] - size2[k-1])
                  / (
                        (size2[k] + size2[k+1])
                      * (size2[k] + size2[k-1])
                    )
                )
              + m3
              * (
                    (size2[k] + size2[k-1])
                  / (
                        (size2[k] + size2[k+1])
                      * (0.5*(size2[k+1] + size2[k-1]) + size2[k])
                    )
                );

            // Write the plane as:
            // sgn(my)*Y = mx*X + mz*Z + alpha*m11*Y = m10*X + m12*Z + alpha

            m1 =
                size1[j-1]*alpha_(i-1,j-1,k)
              + size1[j+1]*alpha_(i-1,j+1,k)
              + size1[j  ]*alpha_(i-1,j,  k);

            m2 =
                size1[j-1]*alpha_(i,j-1,k)
              + size1[j+1]*alpha_(i,j+1,k)
              + size1[j  ]*alpha_(i,j,  k);

            m3 =
                size1[j-1]*alpha_(i+1,j-1,k)
              + size1[j+1]*alpha_(i+1,j+1,k)
              + size1[j  ]*alpha_(i+1,j,  k);

            m[1][0] =
                m1
              * (
                  - (size0[i] + size0[i+1])
                  / (
                        (size0[i] + size0[i-1])
                      * (0.5*(size0[i-1] + size0[i+1]) + size0[i])
                    )
                )
              + m2
              * (
                    2.0*(size0[i+1]-size0[i-1])
                  / (
                        (size0[i] + size0[i+1])
                      * (size0[i] + size0[i-1])
                    )
                )
              + m3
              * (
                    (size0[i] + size0[i-1])
                  / (
                        (size0[i] + size0[i+1])
                      * (0.5*(size0[i+1] + size0[i-1]) + size0[i])
                    )
                );

            m1 =
                size1[j-1]
              * (
                    size0[i  ]*size2[k-1]*alpha_(i,  j-1,k-1)
                  + size0[i  ]*size2[k+1]*alpha_(i,  j-1,k+1)
                  + size0[i-1]*size2[k  ]*alpha_(i-1,j-1,k  )
                  + size0[i+1]*size2[k  ]*alpha_(i+1,j-1,k  )
                  + size0[i  ]*size2[k  ]*alpha_(i,  j-1,k  )
                );

            m2 =
                size1[j+1]
              * (
                    size0[i  ]*size2[k-1]*alpha_(i,  j+1,k-1)
                  + size0[i  ]*size2[k+1]*alpha_(i,  j+1,k+1)
                  + size0[i-1]*size2[k  ]*alpha_(i-1,j+1,k  )
                  + size0[i+1]*size2[k  ]*alpha_(i+1,j+1,k  )
                  + size0[i  ]*size2[k  ]*alpha_(i,  j+1,k  )
                );

            m[1][1] = m1 > m2 ? -1. : 1.;

            m1 =
                size1[j-1]*alpha_(i,j-1,k-1)
              + size1[j+1]*alpha_(i,j+1,k-1)
              + size1[j  ]*alpha_(i,j,  k-1);

            m2 =
                size1[j-1]*alpha_(i,j-1,k)
              + size1[j+1]*alpha_(i,j+1,k)
              + size1[j  ]*alpha_(i,j,  k);

            m3 =
                size1[j-1]*alpha_(i,j-1,k+1)
              + size1[j+1]*alpha_(i,j+1,k+1)
              + size1[j  ]*alpha_(i,j,  k+1);

            m[1][2] =
                m1
              * (
                  - (size2[k] + size2[k+1])
                  / (
                        (size2[k] + size2[k-1])
                      * (0.5*(size2[k-1] + size2[k+1]) + size2[k])
                    )
                )
              + m2
              * (
                    2.0*(size2[k+1] - size2[k-1])
                  / (
                        (size2[k] + size2[k+1])
                      * (size2[k] + size2[k-1])
                    )
                )
              + m3
              * (
                    (size2[k] + size2[k-1])
                  / (
                        (size2[k] + size2[k+1])
                      * (0.5*(size2[k+1] + size2[k-1]) + size2[k])
                    )
                );

            // Write the plane as:
            // sgn(mz)*Z = mx*X + my*Y + alpha*m22*Z = m20*X + m21*Y + alpha

            m1 =
                size2[k-1]*alpha_(i-1,j,k-1)
              + size2[k+1]*alpha_(i-1,j,k+1)
              + size2[k  ]*alpha_(i-1,j,k  );

            m2 =
                size2[k-1]*alpha_(i,j,k-1)
              + size2[k+1]*alpha_(i,j,k+1)
              + size2[k  ]*alpha_(i,j,k  );

            m3 =
                size2[k-1]*alpha_(i+1,j,k-1)
              + size2[k+1]*alpha_(i+1,j,k+1)
              + size2[k  ]*alpha_(i+1,j,k  );

            m[2][0] =
                m1
              * (
                  - (size0[i] + size0[i+1])
                  / (
                        (size0[i] + size0[i-1])
                      * (0.5*(size0[i-1] + size0[i+1]) + size0[i])
                    )
                )
              + m2
              * (
                    2.0*(size0[i+1] - size0[i-1])
                  / (
                        (size0[i] + size0[i+1])
                      * (size0[i] + size0[i-1])
                    )
                )
              + m3
              * (
                    (size0[i] + size0[i-1])
                  / (
                        (size0[i] + size0[i+1])
                      * (0.5*(size0[i+1] + size0[i-1]) + size0[i])
                    )
                );

            m1 =
                size2[k-1]*alpha_(i,j-1,k-1)
              + size2[k+1]*alpha_(i,j-1,k+1)
              + size2[k  ]*alpha_(i,j-1,k  );

            m2 =
                size2[k-1]*alpha_(i,j,k-1)
              + size2[k+1]*alpha_(i,j,k+1)
              + size2[k  ]*alpha_(i,j,k  );

            m3 =
                size2[k-1]*alpha_(i,j+1,k-1)
              + size2[k+1]*alpha_(i,j+1,k+1)
              + size2[k  ]*alpha_(i,j+1,k  );

            m[2][1] =
                m1
              * (
                  - (size1[j] + size1[j+1])
                  / (
                        (size1[j] + size1[j-1])
                      * (0.5*(size1[j-1] + size1[j+1]) + size1[j])
                    )
                )
              + m2
              * (
                    2.0*(size1[j+1] - size1[j-1])
                  / (
                        (size1[j] + size1[j+1])
                      * (size1[j] + size1[j-1])
                    )
                )
              + m3
              * (
                    (size1[j] + size1[j-1])
                  / (
                        (size1[j] + size1[j+1])
                      * (0.5*(size1[j+1] + size1[j-1]) + size1[j])
                    )
                );

            m1 =
                size2[k-1] *
                (
                    size0[i  ]*size1[j-1]*alpha_(i,  j-1,k-1)
                  + size0[i  ]*size1[j+1]*alpha_(i,  j+1,k-1)
                  + size0[i-1]*size1[j  ]*alpha_(i-1,j,  k-1)
                  + size0[i+1]*size1[j  ]*alpha_(i+1,j,  k-1)
                  + size0[i  ]*size1[j  ]*alpha_(i,  j,  k-1)
                );

            m2 =
                size2[k+1]
              * (
                    size0[i  ]*size1[j-1]*alpha_(i,  j-1,k+1)
                  + size0[i  ]*size1[j+1]*alpha_(i,  j+1,k+1)
                  + size0[i-1]*size1[j  ]*alpha_(i-1,j,  k+1)
                  + size0[i+1]*size1[j  ]*alpha_(i+1,j,  k+1)
                  + size0[i  ]*size1[j  ]*alpha_(i,  j,  k+1)
                );

            m[2][2] = m1 > m2 ? -1. : 1.;

            // Normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

            t0 = Foam::mag(m[0][0]) + Foam::mag(m[0][1]) + Foam::mag(m[0][2]);
            m[0][0] /= t0;
            m[0][1] /= t0;
            m[0][2] /= t0;

            t0 = Foam::mag(m[1][0]) + Foam::mag(m[1][1]) + Foam::mag(m[1][2]);
            m[1][0] /= t0;
            m[1][1] /= t0;
            m[1][2] /= t0;

            t0 = Foam::mag(m[2][0]) + Foam::mag(m[2][1]) + Foam::mag(m[2][2]);
            m[2][0] /= t0;
            m[2][1] /= t0;
            m[2][2] /= t0;

            // Choose among the three diagonal values

            t0 = Foam::mag(m[0][0]);
            t1 = Foam::mag(m[1][1]);
            t2 = Foam::mag(m[2][2]);

            label cn = 0;

            if (t1 > t0)
            {
                t0 = t1;
                cn = 1;
            }

            if (t2 > t0)
            {
                cn = 2;
            }

            // Compute gradient averaging the gradient in every vertex of the
            // cell

            double Cx1, Cx2, Cy1, Cy2, Cz1, Cz2;
            m[3][0] = 0.0;
            m[3][1] = 0.0;
            m[3][2] = 0.0;

            for (int il = i - 1; il <= i; il++)
            {
                for (int jl = j - 1; jl <= j; jl++)
                {
                    for (int kl = k - 1; kl <= k; kl++)
                    {
                        Cx1 =
                            (
                                size1[jl  ]*size2[kl  ]*alpha_(il,jl,  kl  )
                              + size1[jl+1]*size2[kl  ]*alpha_(il,jl+1,kl  )
                              + size1[jl  ]*size2[kl+1]*alpha_(il,jl,  kl+1)
                              + size1[jl+1]*size2[kl+1]*alpha_(il,jl+1,kl+1)
                            )
                          / (
                                (size1[jl] + size1[jl+1])
                              * (size2[kl] + size2[kl+1])
                            );

                        Cx2 =
                            (
                                size1[jl  ]*size2[kl  ]*alpha_(il+1,jl,  kl  )
                              + size1[jl+1]*size2[kl  ]*alpha_(il+1,jl+1,kl  )
                              + size1[jl  ]*size2[kl+1]*alpha_(il+1,jl,  kl+1)
                              + size1[jl+1]*size2[kl+1]*alpha_(il+1,jl+1,kl+1)
                            )
                          / (
                                (size1[jl] + size1[jl+1])
                              * (size2[kl] + size2[kl+1])
                            );

                        Cy1 =
                            (
                                size0[il  ]*size2[kl  ]*alpha_(il,  jl,kl  )
                              + size0[il+1]*size2[kl  ]*alpha_(il+1,jl,kl  )
                              + size0[il  ]*size2[kl+1]*alpha_(il,  jl,kl+1)
                              + size0[il+1]*size2[kl+1]*alpha_(il+1,jl,kl+1)
                            )
                          / (
                                (size0[il] + size0[il+1])
                              * (size2[kl] + size2[kl+1])
                            );

                        Cy2 =
                            (
                                size0[il  ]*size2[kl  ]*alpha_(il,  jl+1,kl  )
                              + size0[il+1]*size2[kl  ]*alpha_(il+1,jl+1,kl  )
                              + size0[il  ]*size2[kl+1]*alpha_(il,  jl+1,kl+1)
                              + size0[il+1]*size2[kl+1]*alpha_(il+1,jl+1,kl+1)
                            )
                          / (
                                (size0[il] + size0[il+1])
                              * (size2[kl] + size2[kl+1])
                            );

                        Cz1 =
                            (
                                size0[il  ]*size1[jl  ]*alpha_(il,  jl,  kl)
                              + size0[il+1]*size1[jl  ]*alpha_(il+1,jl,  kl)
                              + size0[il  ]*size1[jl+1]*alpha_(il,  jl+1,kl)
                              + size0[il+1]*size1[jl+1]*alpha_(il+1,jl+1,kl)
                            )
                          / (
                                (size0[il] + size0[il+1])
                              * (size1[jl] + size1[jl+1])
                            );

                        Cz2 =
                            (
                                size0[il  ]*size1[jl  ]*alpha_(il,  jl,  kl+1)
                              + size0[il+1]*size1[jl  ]*alpha_(il+1,jl,  kl+1)
                              + size0[il  ]*size1[jl+1]*alpha_(il,  jl+1,kl+1)
                              + size0[il+1]*size1[jl+1]*alpha_(il+1,jl+1,kl+1)
                            )
                          / (
                                (size0[il] + size0[il+1])
                              * (size1[jl] + size1[jl+1])
                            );

                        m[3][0] += (Cx2 - Cx1)/(16.0*(size0[il] + size0[il+1]));
                        m[3][1] += (Cy2 - Cy1)/(16.0*(size1[jl] + size1[jl+1]));
                        m[3][2] += (Cz2 - Cz1)/(16.0*(size2[kl] + size2[kl+1]));
                    }
                }
            }

            t0 = Foam::mag(m[3][0]) + Foam::mag(m[3][1]) + Foam::mag(m[3][2]);

            if (t0 > 1e-30)
            {
                m[3][0] /= t0;
                m[3][1] /= t0;
                m[3][2] /= t0;

                // Choose between the previous choice and Youngs-CIAM

                t0 = Foam::mag(m[3][0]);
                t1 = Foam::mag(m[3][1]);
                t2 = Foam::mag(m[3][2]);

                if (t1 > t0) t0 = t1;
                if (t2 > t0) t0 = t2;

                if (Foam::mag(m[cn][cn]) < t0)
                    cn = 3;
            }

            // Transform normal to original coordinates

            n(i,j,k) =
                normalised(T.inv() & vector(m[cn][0], m[cn][1], m[cn][2]));
        }
        else
        {
            n(i,j,k) = Zero;
        }
    }

    n.correctBoundaryConditions();
}

}

}

}
