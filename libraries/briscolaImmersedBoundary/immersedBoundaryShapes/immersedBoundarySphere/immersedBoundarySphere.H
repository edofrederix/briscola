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

    // Give the wall distance to the shape
    scalar wallDistance(vector, vector) const;

    // Give the distance normal to the shape
    scalar wallNormalDistance(vector) const;

    // Mirror the given point in the shape surface
    vector mirrorPoint(vector) const;
};

}

}

}

#ifdef NoRepository
    #include "immersedBoundarySphere.C"
#endif

#endif
