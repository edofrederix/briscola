#ifndef level_H
#define level_H

#include "levelPoints.H"
#include "levelBoundaries.H"
#include "levelComms.H"
#include "decomposition.H"

namespace Foam
{

namespace briscola
{

class mesh;

class level
{
    // Reference to the mesh
    const mesh& msh_;

    // Decomposition
    autoPtr<decomposition> decomp_;

    // Level number
    label l_;

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

public:

    // Construct from mesh and parent level (if any)
    level(const mesh&, const level*);

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
        return decomp_;
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

    // Global check if this level is coarsenable
    bool coarsenable() const;

    // Compute this level's coarsened child size
    labelVector coarseN() const;

    // Return parent level if it exists
    const level& parent() const;
};

}

}

#endif
