#ifndef part_H
#define part_H

#include "labelVector.H"
#include "block.H"
#include "face.H"
#include "boundary.H"
#include "partPoints.H"
#include "brickLink.H"

namespace Foam
{

namespace briscola
{

class mesh;

class part
{
    // Level number of this mesh part
    label l_;

    // Size of this mesh part
    labelVector N_;

    // Refinement factors of this level with respect to its parent
    labelVector R_;

    // Mesh points
    partPoints points_;

    // Check if part directions are rectilinear
    labelVector rectilinear_;

    // Check if part directions are uniform
    labelVector uniform_;

    // Bounding box
    faceScalar boundingBox_;

    // Calculate the mesh points
    void calcPoints(const mesh&, const part*);

public:

    // Construct from mesh and part at parent level (if any)
    part(const mesh&, const part*);

    // Copy constructor
    part(const part&);

    // Clone this mesh part
    inline autoPtr<part> clone() const
    {
        return autoPtr<part>(new part(*this));
    }

    // Destructor
    ~part();

    // Level number of this mesh part
    inline label partNum() const { return l_; }

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
    inline const partPoints& points() const { return points_; }

    // Return rectilinear bools
    inline labelVector rectilinear() const { return rectilinear_; }

    // Return uniform bools
    inline labelVector uniform() const { return uniform_; }

    // Return the bounding box
    inline faceScalar boundingBox() const { return boundingBox_; }
};

}

}

#endif
