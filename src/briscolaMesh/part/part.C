#include "part.H"
#include "mesh.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

void part::calcPoints(const mesh& msh, const part* l)
{
    points_.clear();
    points_.setSizeFromCells(N_);

    if (l_ == 0)
    {
        const labelVector start = msh.decomp().myBrickPartStart();

        const brick& b = msh.bricks()[msh.decomp().myBrickNum()];
        const labelVector Nb = b.N();

        // Normalized coordinates of the points in this processor's part of the
        // brick

        scalarField xi(N_.x()+1);
        scalarField eta(N_.y()+1);
        scalarField zeta(N_.z()+1);

        forAll(xi, i)
            xi[i] = scalar(start.x()+i)/Nb.x();

        forAll(eta, j)
            eta[j] = scalar(start.y()+j)/Nb.y();

        forAll(zeta, k)
            zeta[k] = scalar(start.z()+k)/Nb.z();

        const vectorBlock tfi(b.TFI(xi,eta,zeta));

        forAllBlock(tfi, i, j, k)
        {
            points_(i,j,k) = tfi(i,j,k);
        }
    }
    else
    {
        forAllBlock(points_, i, j, k)
        {
            points_(i,j,k) = l->points()(i*R_.x(), j*R_.y(), k*R_.z());
        }
    }

    // Set ghost points

    points_.calcGhostPoints();

    // Remove round-off errors

    points_.clean();
}

part::part(const mesh& msh, const part* l)
:
    l_(l == nullptr ? 0 : l->partNum()+1),
    points_(msh)
{
    if (l == nullptr)
    {
        N_ = msh.decomp().myPartN();
        R_ = zeroXYZ;
    }
    else
    {
        const labelVector P(l->N());

        N_ = msh.coarsen(P);
        R_ = cmptDivide(P,N_);
    }

    calcPoints(msh, l);

    // Check if part directions are rectilinear, which is true if all left,
    // bottom and aft face vectors are aligned

    const scalar tol = 1e-9;

    vector x = zeroXYZ;
    vector y = zeroXYZ;
    vector z = zeroXYZ;

    rectilinear_ = unitXYZ;

    for (int i = 0; i < this->l(); i++)
    for (int j = 0; j < this->m(); j++)
    for (int k = 0; k < this->n(); k++)
    {
        labelVector ijk(i,j,k);

        const vector left =
            0.5
          * (
                (
                    points_(ijk+unitYZ)
                  - points_(ijk)
                )
              ^ (
                    points_(ijk+unitZ)
                  - points_(ijk+unitY)
                )
            );

        const vector bottom =
            0.5
          * (
                (
                    points_(ijk+unitZ)
                  - points_(ijk+unitX)
                )
              ^ (
                    points_(ijk+unitXZ)
                  - points_(ijk)
                )
            );

        const vector aft =
            0.5
          * (
                (
                    points_(ijk+unitXY)
                  - points_(ijk)
                )
              ^ (
                    points_(ijk+unitY)
                  - points_(ijk+unitX)
                )
            );

        if (i == 0 && j == 0 && k == 0)
        {
            x = left;
            y = bottom;
            z = aft;
        }
        else
        {
            if (Foam::mag(x ^ left)/Foam::mag(x) > tol)
            {
                rectilinear_.x() = 0;
            }

            if (Foam::mag(y ^ bottom)/Foam::mag(y) > tol)
            {
                rectilinear_.y() = 0;
            }

            if (Foam::mag(z ^ aft)/Foam::mag(z) > tol)
            {
                rectilinear_.z() = 0;
            }
        }
    }

    // Check if part directions are uniform, which is true if all cell sizes in
    // a direction are the same

    uniform_ = unitXYZ;

    for (int d = 0; d < 3; d++)
    if (rectilinear_[d])
    {
        const vector base = units[d];
        const scalar size = Foam::mag(points_(base) - points_(zeroXYZ));

        for (int i = 1; i < this->N()[d]; i++)
        {
            const scalar next =
                Foam::mag(points_(base*(i+1)) - points_(base*i));

            if (Foam::mag(next-size) > tol)
            {
                uniform_[d] = 0;
                break;
            }
        }
    }
    else
    {
        uniform_[d] = 0;
    }

    // Compute bounding box

    scalarBlock xp(this->N()+unitXYZ);
    scalarBlock yp(this->N()+unitXYZ);
    scalarBlock zp(this->N()+unitXYZ);

    forAllBlock(xp, i, j, k)
    {
        xp(i,j,k) = points_(i,j,k).x();
        yp(i,j,k) = points_(i,j,k).y();
        zp(i,j,k) = points_(i,j,k).z();
    }

    boundingBox_ =
        faceScalar(min(xp), max(xp), min(yp), max(yp), min(zp), max(zp));
}

part::part(const part& l)
:
    l_(l.l_),
    N_(l.N_),
    R_(l.R_),
    points_(l.points_),
    rectilinear_(l.rectilinear_),
    uniform_(l.uniform_),
    boundingBox_(l.boundingBox_)
{}

part::~part()
{}

}

}
