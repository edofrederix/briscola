#ifndef immersedBoundarySlab_H
#define immersedBoundarySlab_H

#include "immersedBoundaryShape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class immersedBoundarySlab
:
    public immersedBoundaryShape
{
    // An infinite slab of given thickness, centered on the plane through
    // center_ with unit normal normal_
    vector center_;
    vector normal_;
    scalar thickness_;

public:

    // Runtime type information
    TypeName("slab");

    // Construct
    immersedBoundarySlab(const dictionary&, bool);

    // Destructor
    virtual ~immersedBoundarySlab();

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
    #include "immersedBoundarySlab.C"
#endif

#endif
