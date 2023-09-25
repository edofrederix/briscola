#ifndef cuboid_H
#define cuboid_H

#include "shape.H"
#include "vector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class cuboid
:
    public shape
{

    // A cuboid is defined by six face centers

    vectorList faceCenters_;

    // The inverted switch can be set to invert cuboid

    bool inverted_;


public:

    // Construct
    cuboid
    (
        vectorList,
        bool
    );

    // Destructor
    virtual ~cuboid();

    bool isInside(vector);

    scalar wallDistance(vector, vector);

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
