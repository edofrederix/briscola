#ifndef immersedBoundaryCylinder_H
#define immersedBoundaryCylinder_H

#include "immersedBoundaryShape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class immersedBoundaryCylinder
:
    public immersedBoundaryShape
{
    // A cylinder is defined by a center line (start, end) and a radius around
    // the center line
    vector start_;
    vector end_;
    scalar radius_;

public:

    // Runtime type information
    TypeName("cylinder");

    // Construct
    immersedBoundaryCylinder(const dictionary&, bool);

    // Destructor
    virtual ~immersedBoundaryCylinder();

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
    #include "immersedBoundaryCylinder.C"
#endif

#endif
