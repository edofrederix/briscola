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
    #include "immersedBoundarySlab.C"
#endif

#endif
