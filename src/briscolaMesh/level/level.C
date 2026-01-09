#include "level.H"
#include "mesh.H"
#include "periodicBoundary.H"

namespace Foam
{

namespace briscola
{

level::level(const mesh& msh, const level* parent)
:
    msh_(msh),
    decomp_(decomposition::New(*this)),
    l_(parent == nullptr ? 0 : parent->levelNum()+1),
    N_
    (
        parent == nullptr
      ? cmptDivide
        (
            msh.bricks()[decomp_->myBrickNum()].N(),
            decomp_->myBrickDecomp()
        )
      : parent->coarseN()
    ),
    R_
    (
        parent == nullptr
      ? zeroXYZ
      : cmptDivide(parent->N(),parent->coarseN())
    ),
    boundaries_(*this),
    comms_(*this),
    points_(*this)
{
    // Check if level directions are rectilinear, which is true if all left,
    // bottom and aft face vectors are aligned

    const scalar tol = 1e-9;

    vector x = zeroXYZ;
    vector y = zeroXYZ;
    vector z = zeroXYZ;

    rectilinear_ = unitXYZ;

    for (int i = 0; i < l(); i++)
    for (int j = 0; j < m(); j++)
    for (int k = 0; k < n(); k++)
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

    // Check if level directions are uniform, which is true if all cell sizes in
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

level::level(const level& lvl)
:
    msh_(lvl.msh_),
    decomp_(lvl.decomp_),
    l_(lvl.l_),
    N_(lvl.N_),
    R_(lvl.R_),
    boundaries_(lvl.boundaries_, *this),
    comms_(lvl.comms_, *this),
    points_(lvl.points_, *this),
    rectilinear_(lvl.rectilinear_),
    uniform_(lvl.uniform_),
    boundingBox_(lvl.boundingBox_)
{}

level::level(const level& lvl, const mesh& msh)
:
    msh_(msh),
    decomp_(lvl.decomp_),
    l_(lvl.l_),
    N_(lvl.N_),
    R_(lvl.R_),
    boundaries_(lvl.boundaries_, *this),
    comms_(lvl.comms_, *this),
    points_(lvl.points_, *this),
    rectilinear_(lvl.rectilinear_),
    uniform_(lvl.uniform_),
    boundingBox_(lvl.boundingBox_)
{}

level::~level()
{}

labelVector level::coarseN() const
{
    if (!msh().topology().structured() || !Pstream::parRun())
    {
        // Keep at least two cells in each direction on unstructured meshes, to
        // avoid heavy distortion

        return labelVector
        (
            N_.x() <= 3 ? N_.x() : N_.x()/2,
            N_.y() <= 3 ? N_.y() : N_.y()/2,
            N_.z() <= 3 ? N_.z() : N_.z()/2
        );
    }
    else
    {
        // Refine until one or three cells in each direction

        labelVector Q
        (
            (N_.x() == 1 || N_.x() == 3) ? N_.x() : N_.x()/2,
            (N_.y() == 1 || N_.y() == 3) ? N_.y() : N_.y()/2,
            (N_.z() == 1 || N_.z() == 3) ? N_.z() : N_.z()/2
        );

        // Avoid having one cell in a periodic direction

        for (int d = 0; d < 3; d++)
        {
            const boundary& b1 = boundaries_.find(faceOffsets[2*d  ]);
            const boundary& b2 = boundaries_.find(faceOffsets[2*d+1]);

            if
            (
                Q[d] < 2
             && b1.castable<periodicBoundary>()
             && b2.castable<periodicBoundary>()
            )
            {
                Q[d] = 2;
            }
        }

        return Q;
    }
}

const level& level::parent() const
{
    #ifdef FULLDEBUG

    if (this->l_ == 0)
        FatalErrorInFunction
            << "Level has no parent" << abort(FatalError);

    #endif

    return msh_[l_-1];
}

}

}
