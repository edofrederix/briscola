#ifndef sphere_H
#define sphere_H

#include "shape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class sphere
:
    public shape
{

    // A sphere is defined by a center point and a radius

    vector center_;

    scalar radius_;

public:

    // Runtime type information
    TypeName("sphere");

    // Construct
    sphere
    (
        const dictionary&,
        bool
    );

    // Destructor
    virtual ~sphere();

    bool isInside(vector);

    scalar wallDistance(vector, vector);

    scalar wallNormalDistance(vector);

    vector mirrorPoint(vector);

};

}

}

}

#ifdef NoRepository
    #include "sphere.C"
#endif

#endif
