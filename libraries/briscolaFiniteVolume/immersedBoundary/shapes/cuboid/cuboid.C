#include "cuboid.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

cuboid::cuboid
(
    vectorList faceCenters,
    bool inverted
)
:
    faceCenters_(faceCenters),
    inverted_(inverted)
{
}

// Destructor

cuboid::~cuboid()
{}

bool cuboid::isInside(vector point)
{
    return false;
}

scalar cuboid::wallDistance(vector c, vector nb)
{
    return -1;
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
