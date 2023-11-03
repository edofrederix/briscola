#ifndef cylinder_H
#define cylinder_H

#include "shape.H"
#include "vector.H"

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

    // The inverted switch can be set to invert cylinder

    bool inverted_;


public:

    // Construct
    cylinder
    (
        vector,
        vector,
        scalar,
        bool
    );

    // Destructor
    virtual ~cylinder();

    bool isInside(vector);

    scalar wallDistance(vector, vector);

    scalar wallNormalDistance(vector);

    vector mirrorPoint(vector);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
