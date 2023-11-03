#ifndef sphere_H
#define sphere_H

#include "shape.H"
#include "vector.H"

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

    // The inverted switch can be set to invert sphere

    bool inverted_;

public:

    // Construct
    sphere
    (
        vector,
        scalar,
        bool
    );

    // Destructor
    virtual ~sphere();

    bool isInside(vector);

    scalar wallDistance(vector, vector);

    scalar wallNormalDistance(vector);

    vector mirrorPoint(vector);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
