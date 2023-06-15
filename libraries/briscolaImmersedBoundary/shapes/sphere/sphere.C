#include "sphere.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

// Constructor

sphere::sphere
(
    vector center,
    scalar radius,
    bool inverted
)
:
    center_(center),
    radius_(radius),
    inverted_(inverted)
{
    if (radius <= 0.0)
    {
        FatalError
            << "Sphere radius has to be a positive value."
            << endl;
        FatalError.exit();
    }
}

// Destructor

sphere::~sphere()
{}

bool sphere::isInside(vector point)
{
    // Check if point is in the plane between the two sphere ends then
    // check f the distance to the sphere axis is smaller than the radius
    if(mag(center_ - point) <= radius_)
    {
        if(!inverted_)
        {
            return true;
        }
    }

    if(!inverted_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
