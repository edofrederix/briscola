#ifndef level_H
#define level_H

#include "levelPoints.H"
#include "levelBoundaries.H"
#include "levelComms.H"
#include "decomposition.H"
#include "labelPair.H"

namespace Foam
{

namespace briscola
{

class mesh;

class level
{
    // Reference to the mesh
    const mesh& msh_;

    // Level number
    label l_;

    // Decomposition
    autoPtr<decomposition> decomp_;

    // Size of this mesh level
    labelVector N_;

    // Refinement factors of this level with respect to its parent
    labelVector R_;

    // Boundaries
    levelBoundaries boundaries_;

    // Parallel communicators
    levelComms comms_;

    // Points
    levelPoints points_;

    // Check if level directions are rectilinear
    labelVector rectilinear_;

    // Check if level directions are uniform
    labelVector uniform_;

    // Bounding box
    faceScalar boundingBox_;

    // Agglomerate mesh size, which can be non-zero in two situations:
    //  1) If we are an agglomerate slave, then the agglomerate mesh size is
    //     that of the agglomerate master. This is needed for the slaves to hold
    //     a copy of the master data.
    //  2) If we are an agglomerate parent and if our child is an agglomerate
    //     master, then the agglomerate mesh size is the sum of all agglomerate
    //     parents. This is needed for the master to hold a copy of the
    //     neighboring data at parent level. This situation is handled by the
    //     child level.
    labelVector Na_;

    // Function to compute this level's N and R from parent pointer
    Pair<labelVector> NandR(const level*);

public:

    // Construct from mesh and parent level (if any)
    level(const mesh&, const level* = nullptr);

    // Copy constructor
    level(const level&);

    // Copy constructor with alternative mesh reference
    level(const level&, const mesh&);

    // Clone this mesh level
    inline autoPtr<level> clone() const
    {
        return autoPtr<level>(new level(*this));
    }

    // Clone this mesh level with alternative mesh reference
    inline autoPtr<level> clone(const mesh& msh) const
    {
        return autoPtr<level>(new level(*this, msh));
    }

    // Destructor
    ~level();

    // Mesh
    inline const mesh& msh() const { return msh_; }

    // Decomposition
    inline const decomposition& decomp() const
    {
        return decomp_();
    }

    // Level number
    inline label levelNum() const { return l_; }

    // Mesh size
    inline labelVector N() const { return N_; }

    // Total number of cells
    inline label size() const { return cmptProduct(N_); }

    // Number of cells along the first direction
    inline label l() const { return N_.x(); }

    // Number of cells along the second direction
    inline label m() const { return N_.y(); }

    // Number of cells along the third direction
    inline label n() const { return N_.z(); }

    // Refinement factors
    inline labelVector R() const { return R_; }

    // Points
    inline const levelPoints& points() const { return points_; }

    // Boundaries
    inline const levelBoundaries& boundaries() const { return boundaries_; }

    // Communicators
    inline const levelComms& comms() const { return comms_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // Return the bounding box
    inline faceScalar boundingBox() const { return boundingBox_; }

    // Mesh size of the agglomerate
    inline labelVector Na() const { return Na_; }

    // Check if this level is coarsenable
    inline bool coarsenable() const
    {
        return N_ != coarsen(N_);
    }

    // Check if level has parent
    inline bool hasParent() const
    {
        return l_ > 0;
    }

    // Return parent level if it exists
    const level& parent() const;

    // Check if level has child
    bool hasChild() const;

    // Return child level if it exists
    const level& child() const;

    // A level is empty if it is not a member of the decomposition
    inline bool empty() const { return !decomp_->member(); }

    // A level is a leaf if it is not empty and either it's on the last level or
    // its child is empty
    inline bool leaf() const
    {
        return !empty() && (!hasChild() || child().empty());
    }

    // Agglomerate parent data offset. When we are an agglomerate slave then we
    // receive a copy of the master's data. The agglomerate parent offset gives
    // the offset to this data which is part of the parent domain.
    inline labelVector aggParentOffset() const
    {
        if (decomp_->aggSlave())
        {
            const labelVector N(cmptDivide(Na_, decomp_->aggProcMap().shape()));
            return cmptMultiply(decomp_->myAggPart(), N);
        }
        else
        {
            return zeroXYZ;
        }
    }
};

}

}

#endif
