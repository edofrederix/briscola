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

MYC::MYC(const vof& vf, const dictionary& dict)
:
    normalScheme(vf, dict),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-3))
{}

MYC::MYC(const MYC& s)
:
    normalScheme(s),
    threshold_(s.threshold_)
{}

MYC::~MYC()
{}

tmp<colocatedVectorField> MYC::operator()()
{
    tmp<colocatedVectorField> tn
    (
        new colocatedVectorField(
            "normal",
            vf_.fvMsh()
        )
    );

    colocatedVectorDirection& n = tn.ref()[0][0];
    n = Zero;

    colocatedScalarField aux = vf_.alpha();
    colocatedScalarDirection alpha = aux[0][0];

    const mesh& mh = vf_.fvMsh().msh();
    const rectilinearMesh& recMh = mh.cast<rectilinearMesh>();
    tensor T = recMh.base();

    PtrList<scalarList> localCellSizes;
    PtrList<scalarList> localPoints;

    const partLevel& part = recMh[0];
    const partLevelPoints& points = part.points();

    localCellSizes.clear();
    localPoints.clear();
    localCellSizes.setSize(3);
    localPoints.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        localCellSizes.set(d, new scalarList(part.N()[d]+2));
        localPoints.set(d, new scalarList(part.N()[d]+3));

        scalarList& localCellSizes_d = localCellSizes[d];
        scalarList& localPoints_d = localPoints[d];

        const labelVector dir = units[d];

        const vector base =
            d == 0 ? T.x() : d == 1 ? T.y() : T.z();

        forAll(localPoints_d, i)
        {
            localPoints_d[i] = points(dir*(i-1)) & base;
        }

        forAll(localCellSizes_d, i)
        {
            localCellSizes_d[i] = localPoints_d[i+1] - localPoints_d[i];
        }
    }

    scalarList xSize = localCellSizes[0];
    scalarList ySize = localCellSizes[1];
    scalarList zSize = localCellSizes[2];

    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > 1e-12) && (alpha(i,j,k) < 1 - 1e-12))
        {
            double m1,m2,m3,m[4][3],t0,t1,t2;
            int cn;
            vector n_aux;

            /* write the plane as: sgn(mx) X =  my Y +  mz Z + alpha
                                        m00 X = m01 Y + m02 Z + alpha */

            m1 = xSize[i] * (ySize[j+1]*zSize[k+1]*alpha(i-1,j,k) + ySize[j+2]*zSize[k+1]*alpha(i-1,j+1,k)
                    + ySize[j]*zSize[k+1]*alpha(i-1,j-1,k) + ySize[j+1]*zSize[k+2]*alpha(i-1,j,k+1) +
                    ySize[j+1]*zSize[k]*alpha(i-1,j,k-1));
            m2 = xSize[i+2] * (ySize[j+1]*zSize[k+1]*alpha(i+1,j,k) + ySize[j+2]*zSize[k+1]*alpha(i+1,j+1,k)
                    + ySize[j]*zSize[k+1]*alpha(i+1,j-1,k) + ySize[j+1]*zSize[k+2]*alpha(i+1,j,k+1) +
                    ySize[j+1]*zSize[k]*alpha(i+1,j,k-1));

            m[0][0] = m1 > m2 ? -1. : 1.;

            m1 = xSize[i]*alpha(i-1,j-1,k) + xSize[i+2]*alpha(i+1,j-1,k) + xSize[i+1]*alpha(i,j-1,k);
            m2 = xSize[i]*alpha(i-1,j,k) + xSize[i+2]*alpha(i+1,j,k) + xSize[i+1]*alpha(i,j,k);
            m3 = xSize[i]*alpha(i-1,j+1,k) + xSize[i+2]*alpha(i+1,j+1,k) + xSize[i+1]*alpha(i,j+1,k);
            m[0][1] = m1 * ((-ySize[j+1]-ySize[j+2])/((ySize[j+1]+ySize[j])*(0.5*(ySize[j]+ySize[j+2])+ySize[j+1])))
                    + m2 * (2*(ySize[j+2]-ySize[j])/((ySize[j+1]+ySize[j+2])*(ySize[j+1]+ySize[j])))
                    + m3 * ((ySize[j+1]+ySize[j])/((ySize[j+1]+ySize[j+2])*(0.5*(ySize[j+2]+ySize[j])+ySize[j+1])));

            m1 = xSize[i]*alpha(i-1,j,k-1) + xSize[i+2]*alpha(i+1,j,k-1) + xSize[i+1]*alpha(i,j,k-1);
            m2 = xSize[i]*alpha(i-1,j,k) + xSize[i+2]*alpha(i+1,j,k) + xSize[i+1]*alpha(i,j,k);
            m3 = xSize[i]*alpha(i-1,j,k+1) + xSize[i+2]*alpha(i+1,j,k+1) + xSize[i+1]*alpha(i,j,k+1);
            m[0][2] = m1 * ((-zSize[k+1]-zSize[k+2])/((zSize[k+1]+zSize[k])*(0.5*(zSize[k]+zSize[k+2])+zSize[k+1])))
                    + m2 * (2*(zSize[k+2]-zSize[k])/((zSize[k+1]+zSize[k+2])*(zSize[k+1]+zSize[k])))
                    + m3 * ((zSize[k+1]+zSize[k])/((zSize[k+1]+zSize[k+2])*(0.5*(zSize[k+2]+zSize[k])+zSize[k+1])));

            /* write the plane as: sgn(my) Y =  mx X +  mz Z + alpha
                                        m11 Y = m10 X + m12 Z + alpha */

            m1 = ySize[j]*alpha(i-1,j-1,k) + ySize[j+2]*alpha(i-1,j+1,k) + ySize[j+1]*alpha(i-1,j,k);
            m2 = ySize[j]*alpha(i,j-1,k) + ySize[j+2]*alpha(i,j+1,k) + ySize[j+1]*alpha(i,j,k);
            m3 = ySize[j]*alpha(i+1,j-1,k) + ySize[j+2]*alpha(i+1,j+1,k) + ySize[j+1]*alpha(i+1,j,k);
            m[1][0] = m1 * ((-xSize[i+1]-xSize[i+2])/((xSize[i+1]+xSize[i])*(0.5*(xSize[i]+xSize[i+2])+xSize[i+1])))
                    + m2 * (2*(xSize[i+2]-xSize[i])/((xSize[i+1]+xSize[i+2])*(xSize[i+1]+xSize[i])))
                    + m3 * ((xSize[i+1]+xSize[i])/((xSize[i+1]+xSize[i+2])*(0.5*(xSize[i+2]+xSize[i])+xSize[i+1])));

            m1 = ySize[j]*(xSize[i+1]*zSize[k]*alpha(i,j-1,k-1) + xSize[i+1]*zSize[k+2]*alpha(i,j-1,k+1)
                    + xSize[i]*zSize[k+1]*alpha(i-1,j-1,k) + xSize[i+2]*zSize[k+1]*alpha(i+1,j-1,k) +
                    xSize[i+1]*zSize[k+1]*alpha(i,j-1,k));
            m2 = ySize[j+2]*(xSize[i+1]*zSize[k]*alpha(i,j+1,k-1) + xSize[i+1]*zSize[k+2]*alpha(i,j+1,k+1)
                    + xSize[i]*zSize[k+1]*alpha(i-1,j+1,k) + xSize[i+2]*zSize[k+1]*alpha(i+1,j+1,k) +
                    xSize[i+1]*zSize[k+1]*alpha(i,j+1,k));
            m[1][1] = m1 > m2 ? -1. : 1.;

            m1 = ySize[j]*alpha(i,j-1,k-1) + ySize[j+2]*alpha(i,j+1,k-1) + ySize[j+1]*alpha(i,j,k-1);
            m2 = ySize[j]*alpha(i,j-1,k) + ySize[j+2]*alpha(i,j+1,k) + ySize[j+1]*alpha(i,j,k);
            m3 = ySize[j]*alpha(i,j-1,k+1) + ySize[j+2]*alpha(i,j+1,k+1) + ySize[j+1]*alpha(i,j,k+1);
            m[1][2] = m1 * ((-zSize[k+1]-zSize[k+2])/((zSize[k+1]+zSize[k])*(0.5*(zSize[k]+zSize[k+2])+zSize[k+1])))
                    + m2 * (2*(zSize[k+2]-zSize[k])/((zSize[k+1]+zSize[k+2])*(zSize[k+1]+zSize[k])))
                    + m3 * ((zSize[k+1]+zSize[k])/((zSize[k+1]+zSize[k+2])*(0.5*(zSize[k+2]+zSize[k])+zSize[k+1])));

            /* write the plane as: sgn(mz) Z =  mx X +  my Y + alpha
                                        m22 Z = m20 X + m21 Y + alpha */

            m1 = zSize[k]*alpha(i-1,j,k-1) + zSize[k+2]*alpha(i-1,j,k+1) + zSize[k+1]*alpha(i-1,j,k);
            m2 = zSize[k]*alpha(i,j,k-1) + zSize[k+2]*alpha(i,j,k+1) + zSize[k+1]*alpha(i,j,k);
            m3 = zSize[k]*alpha(i+1,j,k-1) + zSize[k+2]*alpha(i+1,j,k+1) + zSize[k+1]*alpha(i+1,j,k);
            m[2][0] = m1 * ((-xSize[i+1]-xSize[i+2])/((xSize[i+1]+xSize[i])*(0.5*(xSize[i]+xSize[i+2])+xSize[i+1])))
                    + m2 * (2*(xSize[i+2]-xSize[i])/((xSize[i+1]+xSize[i+2])*(xSize[i+1]+xSize[i])))
                    + m3 * ((xSize[i+1]+xSize[i])/((xSize[i+1]+xSize[i+2])*(0.5*(xSize[i+2]+xSize[i])+xSize[i+1])));

            m1 = zSize[k]*alpha(i,j-1,k-1) + zSize[k+2]*alpha(i,j-1,k+1) + zSize[k+1]*alpha(i,j-1,k);
            m2 = zSize[k]*alpha(i,j,k-1) + zSize[k+2]*alpha(i,j,k+1) + zSize[k+1]*alpha(i,j,k);
            m3 = zSize[k]*alpha(i,j+1,k-1) + zSize[k+2]*alpha(i,j+1,k+1) + zSize[k+1]*alpha(i,j+1,k);
            m[2][1] = m1 * ((-ySize[j+1]-ySize[j+1])/((ySize[j+1]+ySize[j])*(0.5*(ySize[j]+ySize[j+2])+ySize[j+1])))
                    + m2 * (2*(ySize[j+2]-ySize[j])/((ySize[j+1]+ySize[j+2])*(ySize[j+1]+ySize[j])))
                    + m3 * ((ySize[j+1]+ySize[j])/((ySize[j+1]+ySize[j+2])*(0.5*(ySize[j+2]+ySize[j])+ySize[j+1])));

            m1 = zSize[k]*(xSize[i+1]*ySize[j]*alpha(i,j-1,k-1) + xSize[i+1]*ySize[j+2]*alpha(i,j+1,k-1)
                    + xSize[i]*ySize[j+1]*alpha(i-1,j,k-1) + xSize[i+2]*ySize[j+1]*alpha(i+1,j,k-1) +
                    xSize[i+1]*ySize[j+1]*alpha(i,j,k-1));
            m2 = zSize[k+2]*(xSize[i+1]*ySize[j]*alpha(i,j-1,k+1) + xSize[i+1]*ySize[j+2]*alpha(i,j+1,k+1)
                    + xSize[i]*ySize[j+1]*alpha(i-1,j,k+1) + xSize[i+2]*ySize[j+1]*alpha(i+1,j,k+1) +
                    xSize[i+1]*ySize[j+1]*alpha(i,j,k+1));
            m[2][2] = m1 > m2 ? -1. : 1.;

            /* normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1 */

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

            /* choose among the three central scheme */
            t0 = Foam::mag(m[0][0]);
            t1 = Foam::mag(m[1][1]);
            t2 = Foam::mag(m[2][2]);

            cn = 0;
            if (t1 > t0)
            {
                t0 = t1;
                cn = 1;
            }
            if (t2 > t0)
            {
                cn = 2;
            }

            double Cx1, Cx2, Cy1, Cy2, Cz1, Cz2;
            m[3][0] = 0;
            m[3][1] = 0;
            m[3][2] = 0;

            for (int il = i - 1; il <= i; il++)
            {
                for (int jl = j - 1; jl <= j; jl++)
                {
                    for (int kl = k - 1; kl <= k; kl++)
                    {
                        Cx1 = (ySize[jl+1]*zSize[kl+1]*alpha(il,jl,kl) + ySize[jl+2]*zSize[kl+1]*alpha(il,jl+1,kl)
                                + ySize[jl+1]*zSize[kl+2]*alpha(il,jl,kl+1) + ySize[jl+2]*zSize[kl+2]*alpha(il,jl+1,kl+1))
                                /((ySize[jl+1] + ySize[jl + 2])*(zSize[kl+1] + zSize[kl + 2]));
                        Cx2 = (ySize[jl+1]*zSize[kl+1]*alpha(il+1,jl,kl) + ySize[jl+2]*zSize[kl+1]*alpha(il+1,jl+1,kl)
                                + ySize[jl+1]*zSize[kl+2]*alpha(il+1,jl,kl+1) + ySize[jl+2]*zSize[kl+2]*alpha(il+1,jl+1,kl+1))
                                /((ySize[jl+1] + ySize[jl + 2])*(zSize[kl+1] + zSize[kl + 2]));
                        Cy1 = (xSize[il+1]*zSize[kl+1]*alpha(il,jl,kl) + xSize[il+2]*zSize[kl+1]*alpha(il+1,jl,kl)
                                + xSize[il+1]*zSize[kl+2]*alpha(il,jl,kl+1) + xSize[il+2]*zSize[kl+2]*alpha(il+1,jl,kl+1))
                                /((xSize[il+1] + xSize[il + 2])*(zSize[kl+1] + zSize[kl + 2]));
                        Cy2 = (xSize[il+1]*zSize[kl+1]*alpha(il,jl+1,kl) + xSize[il+2]*zSize[kl+1]*alpha(il+1,jl+1,kl)
                                + xSize[il+1]*zSize[kl+2]*alpha(il,jl+1,kl+1) + xSize[il+2]*zSize[kl+2]*alpha(il+1,jl+1,kl+1))
                                /((xSize[il+1] + xSize[il + 2])*(zSize[kl+1] + zSize[kl + 2]));
                        Cz1 = (xSize[il+1]*ySize[jl+1]*alpha(il,jl,kl) + xSize[il+2]*ySize[jl+1]*alpha(il+1,jl,kl)
                                + xSize[il+1]*ySize[jl+2]*alpha(il,jl+1,kl) + xSize[il+2]*ySize[jl+2]*alpha(il+1,jl+1,kl))
                                /((xSize[il+1] + xSize[il + 2])*(ySize[jl+1] + ySize[jl + 2]));
                        Cz2 = (xSize[il+1]*ySize[jl+1]*alpha(il,jl,kl+1) + xSize[il+2]*ySize[jl+1]*alpha(il+1,jl,kl+1)
                                + xSize[il+1]*ySize[jl+2]*alpha(il,jl+1,kl+1) + xSize[il+2]*ySize[jl+2]*alpha(il+1,jl+1,kl+1))
                                /((xSize[il+1] + xSize[il + 2])*(ySize[jl+1] + ySize[jl + 2]));

                        m[3][0] += (Cx2 - Cx1)/(16 * (xSize[il+1] + xSize[il+2]));
                        m[3][1] += (Cy2 - Cy1)/(16 * (ySize[jl+1] + ySize[jl+2]));
                        m[3][2] += (Cz2 - Cz1)/(16 * (zSize[kl+1] + zSize[kl+2]));

                    }
                }
            }

            t0 = Foam::mag(m[3][0]) + Foam::mag(m[3][1]) + Foam::mag(m[3][2]);
            if (t0 > 1e-30)
            {
                m[3][0] /= t0;
                m[3][1] /= t0;
                m[3][2] /= t0;

                /* choose between the previous choice and Youngs-CIAM */
                t0 = Foam::mag(m[3][0]);
                t1 = Foam::mag(m[3][1]);
                t2 = Foam::mag(m[3][2]);
                if (t1 > t0)
                    t0 = t1;
                if (t2 > t0)
                    t0 = t2;

                if (Foam::mag(m[cn][cn]) < t0)
                    cn = 3;
            }

            n_aux[0] = m[cn][0];
            n_aux[1] = m[cn][1];
            n_aux[2] = m[cn][2];

            n_aux = T.inv() & n_aux;
            n(i,j,k) = n_aux / Foam::mag(n_aux);

        }
    }

    return tn;
}

}

}

}
