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
    threshold_(vf.threshold())
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

    /*

    +++++++++++++++++++++++++++++++++++++++++++++++

    Reconstruct the normal interface usign the Mixed
    Youngs method in Aulisa (2007).

    The method has been generalized for non uniform
    meshes using three points lagrange derivative
    aproximations and for skwewed meshes by transforming
    the coorinate system into a perpendicular one.

    +++++++++++++++++++++++++++++++++++++++++++++++

    */

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

    PtrList<PartialList<scalar>> localCellSizes = recMh.localCellSizes();

    PartialList<scalar> xSize = localCellSizes[0];
    PartialList<scalar> ySize = localCellSizes[1];
    PartialList<scalar> zSize = localCellSizes[2];

    forAllCells(n, i, j, k)
    {
        if ((alpha(i,j,k) > threshold_) && (alpha(i,j,k) < 1 - threshold_))
        {
            double m1,m2,m3,m[4][3],t0,t1,t2;
            int cn;
            vector n_aux;

            /* write the plane as: sgn(mx) X =  my Y +  mz Z + alpha
                                        m00 X = m01 Y + m02 Z + alpha */

            m1 = xSize[i-1] * (ySize[j]*zSize[k]*alpha(i-1,j,k) + ySize[j+1]*zSize[k]*alpha(i-1,j+1,k)
                    + ySize[j-1]*zSize[k]*alpha(i-1,j-1,k) + ySize[j]*zSize[k+1]*alpha(i-1,j,k+1) +
                    ySize[j]*zSize[k-1]*alpha(i-1,j,k-1));
            m2 = xSize[i+1] * (ySize[j]*zSize[k]*alpha(i+1,j,k) + ySize[j+1]*zSize[k]*alpha(i+1,j+1,k)
                    + ySize[j-1]*zSize[k]*alpha(i+1,j-1,k) + ySize[j]*zSize[k+1]*alpha(i+1,j,k+1) +
                    ySize[j]*zSize[k-1]*alpha(i+1,j,k-1));

            m[0][0] = m1 > m2 ? -1. : 1.;

            m1 = xSize[i-1]*alpha(i-1,j-1,k) + xSize[i+1]*alpha(i+1,j-1,k) + xSize[i]*alpha(i,j-1,k);
            m2 = xSize[i-1]*alpha(i-1,j,k) + xSize[i+1]*alpha(i+1,j,k) + xSize[i]*alpha(i,j,k);
            m3 = xSize[i-1]*alpha(i-1,j+1,k) + xSize[i+1]*alpha(i+1,j+1,k) + xSize[i]*alpha(i,j+1,k);
            m[0][1] = m1 * ((-ySize[j]-ySize[j+1])/((ySize[j]+ySize[j-1])*(0.5*(ySize[j-1]+ySize[j+1])+ySize[j])))
                    + m2 * (2*(ySize[j+1]-ySize[j-1])/((ySize[j]+ySize[j+1])*(ySize[j]+ySize[j-1])))
                    + m3 * ((ySize[j]+ySize[j-1])/((ySize[j]+ySize[j+1])*(0.5*(ySize[j+1]+ySize[j-1])+ySize[j])));

            m1 = xSize[i-1]*alpha(i-1,j,k-1) + xSize[i+1]*alpha(i+1,j,k-1) + xSize[i]*alpha(i,j,k-1);
            m2 = xSize[i-1]*alpha(i-1,j,k) + xSize[i+1]*alpha(i+1,j,k) + xSize[i]*alpha(i,j,k);
            m3 = xSize[i-1]*alpha(i-1,j,k+1) + xSize[i+1]*alpha(i+1,j,k+1) + xSize[i]*alpha(i,j,k+1);
            m[0][2] = m1 * ((-zSize[k]-zSize[k+1])/((zSize[k]+zSize[k-1])*(0.5*(zSize[k-1]+zSize[k+1])+zSize[k])))
                    + m2 * (2*(zSize[k+1]-zSize[k-1])/((zSize[k]+zSize[k+1])*(zSize[k]+zSize[k-1])))
                    + m3 * ((zSize[k]+zSize[k-1])/((zSize[k]+zSize[k+1])*(0.5*(zSize[k+1]+zSize[k-1])+zSize[k])));

            /* write the plane as: sgn(my) Y =  mx X +  mz Z + alpha
                                        m11 Y = m10 X + m12 Z + alpha */

            m1 = ySize[j-1]*alpha(i-1,j-1,k) + ySize[j+1]*alpha(i-1,j+1,k) + ySize[j]*alpha(i-1,j,k);
            m2 = ySize[j-1]*alpha(i,j-1,k) + ySize[j+1]*alpha(i,j+1,k) + ySize[j]*alpha(i,j,k);
            m3 = ySize[j-1]*alpha(i+1,j-1,k) + ySize[j+1]*alpha(i+1,j+1,k) + ySize[j]*alpha(i+1,j,k);
            m[1][0] = m1 * ((-xSize[i]-xSize[i+1])/((xSize[i]+xSize[i-1])*(0.5*(xSize[i-1]+xSize[i+1])+xSize[i])))
                    + m2 * (2*(xSize[i+1]-xSize[i-1])/((xSize[i]+xSize[i+1])*(xSize[i]+xSize[i-1])))
                    + m3 * ((xSize[i]+xSize[i-1])/((xSize[i]+xSize[i+1])*(0.5*(xSize[i+1]+xSize[i-1])+xSize[i])));

            m1 = ySize[j-1]*(xSize[i]*zSize[k-1]*alpha(i,j-1,k-1) + xSize[i]*zSize[k+1]*alpha(i,j-1,k+1)
                    + xSize[i-1]*zSize[k]*alpha(i-1,j-1,k) + xSize[i+1]*zSize[k]*alpha(i+1,j-1,k) +
                    xSize[i]*zSize[k]*alpha(i,j-1,k));
            m2 = ySize[j+1]*(xSize[i]*zSize[k-1]*alpha(i,j+1,k-1) + xSize[i]*zSize[k+1]*alpha(i,j+1,k+1)
                    + xSize[i-1]*zSize[k]*alpha(i-1,j+1,k) + xSize[i+1]*zSize[k]*alpha(i+1,j+1,k) +
                    xSize[i]*zSize[k]*alpha(i,j+1,k));
            m[1][1] = m1 > m2 ? -1. : 1.;

            m1 = ySize[j-1]*alpha(i,j-1,k-1) + ySize[j+1]*alpha(i,j+1,k-1) + ySize[j]*alpha(i,j,k-1);
            m2 = ySize[j-1]*alpha(i,j-1,k) + ySize[j+1]*alpha(i,j+1,k) + ySize[j]*alpha(i,j,k);
            m3 = ySize[j-1]*alpha(i,j-1,k+1) + ySize[j+1]*alpha(i,j+1,k+1) + ySize[j]*alpha(i,j,k+1);
            m[1][2] = m1 * ((-zSize[k]-zSize[k+1])/((zSize[k]+zSize[k-1])*(0.5*(zSize[k-1]+zSize[k+1])+zSize[k])))
                    + m2 * (2*(zSize[k+1]-zSize[k-1])/((zSize[k]+zSize[k+1])*(zSize[k]+zSize[k-1])))
                    + m3 * ((zSize[k]+zSize[k-1])/((zSize[k]+zSize[k+1])*(0.5*(zSize[k+1]+zSize[k-1])+zSize[k])));

            /* write the plane as: sgn(mz) Z =  mx X +  my Y + alpha
                                        m22 Z = m20 X + m21 Y + alpha */

            m1 = zSize[k-1]*alpha(i-1,j,k-1) + zSize[k+1]*alpha(i-1,j,k+1) + zSize[k]*alpha(i-1,j,k);
            m2 = zSize[k-1]*alpha(i,j,k-1) + zSize[k+1]*alpha(i,j,k+1) + zSize[k]*alpha(i,j,k);
            m3 = zSize[k-1]*alpha(i+1,j,k-1) + zSize[k+1]*alpha(i+1,j,k+1) + zSize[k]*alpha(i+1,j,k);
            m[2][0] = m1 * ((-xSize[i]-xSize[i+1])/((xSize[i]+xSize[i-1])*(0.5*(xSize[i-1]+xSize[i+1])+xSize[i])))
                    + m2 * (2*(xSize[i+1]-xSize[i-1])/((xSize[i]+xSize[i+1])*(xSize[i]+xSize[i-1])))
                    + m3 * ((xSize[i]+xSize[i-1])/((xSize[i]+xSize[i+1])*(0.5*(xSize[i+1]+xSize[i-1])+xSize[i])));

            m1 = zSize[k-1]*alpha(i,j-1,k-1) + zSize[k+1]*alpha(i,j-1,k+1) + zSize[k]*alpha(i,j-1,k);
            m2 = zSize[k-1]*alpha(i,j,k-1) + zSize[k+1]*alpha(i,j,k+1) + zSize[k]*alpha(i,j,k);
            m3 = zSize[k-1]*alpha(i,j+1,k-1) + zSize[k+1]*alpha(i,j+1,k+1) + zSize[k]*alpha(i,j+1,k);
            m[2][1] = m1 * ((-ySize[j]-ySize[j+1])/((ySize[j]+ySize[j-1])*(0.5*(ySize[j-1]+ySize[j+1])+ySize[j])))
                    + m2 * (2*(ySize[j+1]-ySize[j-1])/((ySize[j]+ySize[j+1])*(ySize[j]+ySize[j-1])))
                    + m3 * ((ySize[j]+ySize[j-1])/((ySize[j]+ySize[j+1])*(0.5*(ySize[j+1]+ySize[j-1])+ySize[j])));

            m1 = zSize[k-1]*(xSize[i]*ySize[j-1]*alpha(i,j-1,k-1) + xSize[i]*ySize[j+1]*alpha(i,j+1,k-1)
                    + xSize[i-1]*ySize[j]*alpha(i-1,j,k-1) + xSize[i+1]*ySize[j]*alpha(i+1,j,k-1) +
                    xSize[i]*ySize[j]*alpha(i,j,k-1));
            m2 = zSize[k+1]*(xSize[i]*ySize[j-1]*alpha(i,j-1,k+1) + xSize[i]*ySize[j+1]*alpha(i,j+1,k+1)
                    + xSize[i-1]*ySize[j]*alpha(i-1,j,k+1) + xSize[i+1]*ySize[j]*alpha(i+1,j,k+1) +
                    xSize[i]*ySize[j]*alpha(i,j,k+1));
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

            /* compute gradient averaging the gradient in every vertex of the cell*/

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
                        Cx1 = (ySize[jl]*zSize[kl]*alpha(il,jl,kl) + ySize[jl+1]*zSize[kl]*alpha(il,jl+1,kl)
                                + ySize[jl]*zSize[kl+1]*alpha(il,jl,kl+1) + ySize[jl+1]*zSize[kl+1]*alpha(il,jl+1,kl+1))
                                /((ySize[jl] + ySize[jl+1])*(zSize[kl] + zSize[kl+1]));
                        Cx2 = (ySize[jl]*zSize[kl]*alpha(il+1,jl,kl) + ySize[jl+1]*zSize[kl]*alpha(il+1,jl+1,kl)
                                + ySize[jl]*zSize[kl+1]*alpha(il+1,jl,kl+1) + ySize[jl+1]*zSize[kl+1]*alpha(il+1,jl+1,kl+1))
                                /((ySize[jl] + ySize[jl+1])*(zSize[kl] + zSize[kl+1]));
                        Cy1 = (xSize[il]*zSize[kl]*alpha(il,jl,kl) + xSize[il+1]*zSize[kl]*alpha(il+1,jl,kl)
                                + xSize[il]*zSize[kl+1]*alpha(il,jl,kl+1) + xSize[il+1]*zSize[kl+1]*alpha(il+1,jl,kl+1))
                                /((xSize[il] + xSize[il+1])*(zSize[kl] + zSize[kl+1]));
                        Cy2 = (xSize[il]*zSize[kl]*alpha(il,jl+1,kl) + xSize[il+1]*zSize[kl]*alpha(il+1,jl+1,kl)
                                + xSize[il]*zSize[kl+1]*alpha(il,jl+1,kl+1) + xSize[il+1]*zSize[kl+1]*alpha(il+1,jl+1,kl+1))
                                /((xSize[il] + xSize[il+1])*(zSize[kl] + zSize[kl+1]));
                        Cz1 = (xSize[il]*ySize[jl]*alpha(il,jl,kl) + xSize[il+1]*ySize[jl]*alpha(il+1,jl,kl)
                                + xSize[il]*ySize[jl+1]*alpha(il,jl+1,kl) + xSize[il+1]*ySize[jl+1]*alpha(il+1,jl+1,kl))
                                /((xSize[il] + xSize[il+1])*(ySize[jl] + ySize[jl+1]));
                        Cz2 = (xSize[il]*ySize[jl]*alpha(il,jl,kl+1) + xSize[il+1]*ySize[jl]*alpha(il+1,jl,kl+1)
                                + xSize[il]*ySize[jl+1]*alpha(il,jl+1,kl+1) + xSize[il+1]*ySize[jl+1]*alpha(il+1,jl+1,kl+1))
                                /((xSize[il] + xSize[il+1])*(ySize[jl] + ySize[jl+1]));

                        m[3][0] += (Cx2 - Cx1)/(16 * (xSize[il] + xSize[il+1]));
                        m[3][1] += (Cy2 - Cy1)/(16 * (ySize[jl] + ySize[jl+1]));
                        m[3][2] += (Cz2 - Cz1)/(16 * (zSize[kl] + zSize[kl+1]));

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

            /* set normal */
            n_aux[0] = m[cn][0];
            n_aux[1] = m[cn][1];
            n_aux[2] = m[cn][2];

            /* transform normal to orginal coordinates*/
            n_aux = T.inv() & n_aux;
            n(i,j,k) = n_aux / Foam::mag(n_aux);

        }
    }

    return tn;
}

}

}

}
