#ifndef immersedBoundarySphere_H
#define immersedBoundarySphere_H

#include "immersedBoundaryShape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class immersedBoundarySphere
:
    public immersedBoundaryShape
{
    // A sphere is defined by a center point and a radius
    vector center_;
    scalar radius_;

public:

    // Runtime type information
    TypeName("sphere");

    // Construct
    immersedBoundarySphere(const dictionary&, bool);

    // Destructor
    virtual ~immersedBoundarySphere();

    // Return true of the vector is inside the shape
    bool isInside(vector) const;

    // Return dimensionless wall distance given two points. The first point must
    // be in the fluid and the second point in the IB.
    virtual scalar wallDistance(vector, vector) const;

    // Return the signed distance from a point to the immersed boundary surface
    // in the wall-normal direction. Negative means the point is in the fluid.
    virtual scalar wallNormalDistance(vector) const;

    // Returns the coordinates of the mirror point in the wall-normal direction
    virtual vector mirrorPoint(vector) const;
};

}

}

}

#ifdef NoRepository
    #include "immersedBoundarySphere.C"
#endif

#endif
