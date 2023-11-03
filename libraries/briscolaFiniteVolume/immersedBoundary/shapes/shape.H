#ifndef shape_H
#define shape_H

#include "vector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class shape
{

public:

    // Empty constructor
    shape();

    // Destructor
    virtual ~shape();

    // Check if point is inside shape
    virtual bool isInside(vector) = 0;

    // Returns dimensionless wall distance given two points
    // The first point is a fluid cell and the second point is an IB cell
    virtual scalar wallDistance(vector, vector) = 0;

    // Returns the distance from a ghost cell to the immersed boundary in
    // the wall-normal direction
    virtual scalar wallNormalDistance(vector) = 0;

    // Returns the coordinates of the mirror point of a ghost cell
    // in the wall-normal direction of a shape
    virtual vector mirrorPoint(vector) = 0;

};

} // end namespace fv

} // end namespace briscola

} // end namespace Foam

#endif
