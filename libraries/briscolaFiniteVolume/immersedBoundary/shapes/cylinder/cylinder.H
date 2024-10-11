#ifndef cylinder_H
#define cylinder_H

#include "shape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class cylinder
:
    public shape
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
    cylinder(const dictionary&, bool);

    // Destructor
    virtual ~cylinder();

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
    #include "cylinder.C"
#endif

#endif
