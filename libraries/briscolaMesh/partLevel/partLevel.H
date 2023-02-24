#ifndef partLevel_H
#define partLevel_H

#include "labelVector.H"
#include "block.H"
#include "face.H"
#include "partPatch.H"
#include "partLevelPoints.H"
#include "brickLink.H"

namespace Foam
{

namespace briscola
{

class mesh;

class partLevel
{
    // Reference to mesh object
    const mesh& msh_;

    // Pointer to the parent mesh part level
    const partLevel* parentPtr_;

    // Part level number
    label l_;

    // Mesh size of this part level
    labelVector N_;

    // Refinement factors of this part level with respect to its parent
    labelVector R_;

    // Mesh points
    partLevelPoints points_;

    // Check if part directions are rectilinear
    labelVector rectilinear_;

    // Check if part directions are uniform
    labelVector uniform_;

    // Calculate the mesh points
    void calcPoints();

    // Calculate ghost points of the mesh
    void calcGhostPoints();

public:

    // Construct from mesh and parent part level (if any)
    partLevel(const mesh&, const partLevel*);

    // Destructor
    ~partLevel();

    // Reference to the mesh
    inline const mesh& msh() const { return msh_; }

    // Reference to the parent mesh part level (unchecked!)
    inline const partLevel& parent() const { return *parentPtr_; }

    // Pointer to the parent mesh part level
    inline const partLevel* parentPtr() { return parentPtr_; }

    // Part level number
    inline label partLevelNum() const { return l_; }

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
    inline const partLevelPoints& points() const { return points_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // For rectlinear directions, return the cell sizes. For non-rectilinear
    // directions, return a zero list.
    scalarList rectilinearCellSizes(const label) const;

    // Return the start cell of the internal cell boundary at given offset
    inline labelVector boundaryStart(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? 0 : (bo.x()+1)/2*(l()-1),
            bo.y() == 0 ? 0 : (bo.y()+1)/2*(m()-1),
            bo.z() == 0 ? 0 : (bo.z()+1)/2*(n()-1)
        );
    }

    // Return the shape of the internal cell boundary at given offset
    inline labelVector boundaryN(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? l() : 1,
            bo.y() == 0 ? m() : 1,
            bo.z() == 0 ? n() : 1
        );
    }

    // Return the end of the internal cell boundary at given offset
    inline labelVector boundaryEnd(const labelVector& bo) const
    {
        return boundaryStart(bo) + boundaryN(bo);
    }

    // Return the start point of the point boundary at given offset
    inline labelVector pointBoundaryStart(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? 0 : (bo.x()+1)/2*l(),
            bo.y() == 0 ? 0 : (bo.y()+1)/2*m(),
            bo.z() == 0 ? 0 : (bo.z()+1)/2*n()
        );
    }

    // Return the shape of the point boundary at given offset
    inline labelVector pointBoundaryN(const labelVector& bo) const
    {
        return labelVector
        (
            bo.x() == 0 ? l()+1 : 1,
            bo.y() == 0 ? m()+1 : 1,
            bo.z() == 0 ? n()+1 : 1
        );
    }

    // Return the end of the point boundary at given offset
    inline labelVector pointBoundaryEnd(const labelVector& bo) const
    {
        return pointBoundaryStart(bo) + pointBoundaryN(bo);
    }
};

}

}

#endif
