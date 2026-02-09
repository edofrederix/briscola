#include "level.H"
#include "mesh.H"
#include "periodicBoundary.H"

namespace Foam
{

namespace briscola
{

Pair<labelVector> level::NandR(const level* parent)
{
    labelVector N, R;

    if (parent == nullptr)
    {
        // Top level, just divide the brick's N by the brick decomposition

        N =
            cmptDivide
            (
                decomp_->lvl().msh().bricks()[decomp_->myBrickNum()].N(),
                decomp_->myBrickDecomp()
            );

        R = zeroXYZ;
    }
    else
    {
        if (decomp_->agglomerated())
        {
            const labelVector NR = cmptMultiply(parent->N(), decomp_->R());

            N = coarsen(NR);
            R = cmptDivide(NR, N);

            if (decomp_->aggSlave())
                N = Zero;
        }
        else if (decomp_->member())
        {
            N = coarsen(parent->N());
            R = cmptDivide(parent->N(), N);
        }
        else
        {
            N = Zero;
            R = Zero;
        }
    }

    return Pair<labelVector>(N,R);
}

level::level(const mesh& msh, const level* parentPtr)
:
    msh_(msh),
    l_(parentPtr == nullptr ? 0 : parentPtr->levelNum()+1),
    decomp_(decomposition::New(*this)),
    N_(this->NandR(parentPtr).first()),
    R_(this->NandR(parentPtr).second()),
    boundaries_(*this),
    comms_(*this),
    points_(*this),
    rectilinear_(Zero),
    uniform_(Zero),
    boundingBox_(Zero),
    Na_(Zero)
{
    // The agglomerate mesh size is non-zero if we are an agglomerate slave.
    // Receive from master.

    if (decomp_->aggSlave())
    {
        // Receive from master

        IPstream fromMaster
        (
            Pstream::commsTypes::blocking,
            0,
            0,
            Pstream::msgType(),
            comms_.agg()
        );

        fromMaster >> Na_;
    }
    else if (decomp_->aggMaster())
    {
        // Send to slaves

        forAllBlockLinear(decomp_->aggProcMap(), proc)
        {
            if (proc > 0)
            {
                OPstream toSlave
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    Pstream::msgType(),
                    comms_.agg()
                );

                toSlave << N_;
            }
        }
    }

    // If we are an agglomerate member, then our parent is an agglomerate
    // parent. In that case the master process has an agglomerate parent mesh
    // size which is the sum of all agglomerate parent mesh sizes. We modify the
    // parent. This is done by the child because upon construction of the parent
    // there's no child information available yet.

    if (decomp_->aggSlave())
    {
        // Send to master

        OPstream toMaster
        (
            Pstream::commsTypes::blocking,
            0,
            0,
            Pstream::msgType(),
            comms_.agg()
        );

        toMaster << parent().N();
    }
    else if (decomp_->aggMaster())
    {
        // Receive from slaves

        level& parent = const_cast<level&>(*parentPtr);

        labelVectorBlock Ns(decomp_->aggProcMap().shape());
        Ns(0,0,0) = parent.N();

        label proc = 0;
        forAllBlock(Ns, i, j, k)
        {
            if (proc != 0)
            {
                IPstream fromSlave
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    Pstream::msgType(),
                    comms_.agg()
                );

                labelVector Na;
                fromSlave >> Na;

                Ns(i,j,k) = Na;
            }

            proc++;
        }

        // Sum sizes and store in parent

        parent.Na_ = Zero;

        for (int d = 0; d < 3; d++)
            for (int i = 0; i < Ns.shape()[d]; i++)
                parent.Na_[d] += Ns(units[d]*i)[d];

        if (cmptProduct(parent.N_) >= cmptProduct(parent.Na_))
            FatalErrorInFunction
                << "Invalid agglomerate size for parent. "
                << "Parent size = " << parent.N_ << ", "
                << "Agglomerate size = " << parent.Na_ << endl
                << abort(FatalError);
    }

    // No need to further initialize empty levels

    if (empty())
        return;

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
            #ifdef FULLDEBUG

            if
            (
                Foam::mag(x) < VSMALL
             || Foam::mag(y) < VSMALL
             || Foam::mag(z) < VSMALL
            )
            {
                FatalErrorInFunction
                    << "Degenerate cell detected at " << ijk
                    << abort(FatalError);
            }

            #endif

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
    l_(lvl.l_),
    decomp_(decomposition::New(*this)),
    N_(lvl.N_),
    R_(lvl.R_),
    boundaries_(lvl.boundaries_, *this),
    comms_(lvl.comms_, *this),
    points_(lvl.points_, *this),
    rectilinear_(lvl.rectilinear_),
    uniform_(lvl.uniform_),
    boundingBox_(lvl.boundingBox_),
    Na_(lvl.Na_)
{}

level::level(const level& lvl, const mesh& msh)
:
    msh_(msh),
    l_(lvl.l_),
    decomp_(decomposition::New(*this)),
    N_(lvl.N_),
    R_(lvl.R_),
    boundaries_(lvl.boundaries_, *this),
    comms_(lvl.comms_, *this),
    points_(lvl.points_, *this),
    rectilinear_(lvl.rectilinear_),
    uniform_(lvl.uniform_),
    boundingBox_(lvl.boundingBox_),
    Na_(lvl.Na_)
{}

level::~level()
{}

const level& level::parent() const
{
    #ifdef FULLDEBUG

    if (!hasParent())
        FatalErrorInFunction
            << "Level has no parent" << abort(FatalError);

    #endif

    return msh_[l_-1];
}

bool level::hasChild() const
{
    return l_ < msh_.size() - 1;
}

const level& level::child() const
{
    #ifdef FULLDEBUG
    if (l_ == msh_.size()-1)
        FatalErrorInFunction
            << "Level has no child" << abort(FatalError);
    #endif

    return msh_[l_+1];
}

}

}
