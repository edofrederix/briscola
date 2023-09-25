#include "sphere.H"

namespace Foam
{

namespace briscola
{

namespace fv
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
    // Check if distance from point to sphere center
    // is smaller than the sphere's radius
    if(mag(center_ - point) <= radius_)
    {
        if(!inverted_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else if(!inverted_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

scalar sphere::wallDistance(vector c, vector nb)
{
    // Return -1 if the center point is not a fluid point
    // or if the neighboring point is not inside the sphere
    if (this->isInside(c))
    {
        return -1;
    }

    if (!this->isInside(nb))
    {
        return -1;
    }

    // Normalized direction vector of the line
    vector D = (nb-c)/mag(nb-c);

    // Vector from origin of the line to center of the sphere
    vector L = center_-c;

    scalar tc = L & D;

    if (tc < 0)
    {
        return -1;
    }

    scalar d = sqrt(magSqr(L)-sqr(tc));

    if (d > radius_)
    {
        return -1;
    }

    scalar t1c = sqrt(sqr(radius_) - sqr(d));

    return (tc - t1c);
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
