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

    // A cylinder is defined by a center line (start, end)
    // and a radius around the center line

    vector start_;

    vector end_;

    scalar radius_;

public:

    // Runtime type information
    TypeName("cylinder");

    // Construct
    cylinder
    (
        const dictionary&,
        bool
    );

    // Destructor
    virtual ~cylinder();

    bool isInside(vector);

    scalar wallDistance(vector, vector);

    scalar wallNormalDistance(vector);

    vector mirrorPoint(vector);

};

}

}

}

#ifdef NoRepository
    #include "cylinder.C"
#endif

#endif
