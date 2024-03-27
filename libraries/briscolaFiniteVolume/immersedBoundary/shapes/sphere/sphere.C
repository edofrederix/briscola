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
    const dictionary& dict,
    bool inverted
)
:
    shape(dict,inverted),
    center_(vector(dict.lookup("center"))),
    radius_(readScalar(dict.lookup("radius")))
{
    if (radius_ <= 0.0)
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
        if(!this->inverted_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else if(!this->inverted_)
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

scalar sphere::wallNormalDistance(vector gc)
{
    // Return -1 if the point is outside of the sphere
    if (!this->isInside(gc))
    {
        return -1;
    }

    // Return radius minus distance from center to ghost cell
    return (inverted_ ? mag(gc-center_) - radius_ : radius_ - mag(gc-center_));
}

vector sphere::mirrorPoint(vector gc)
{
    scalar dist = this->wallNormalDistance(gc);

    if (dist < 0)
    {
        // This could be a problem if gc is exactly on the IB
        return gc;
    }

    // Wall-normal unit vector
    vector n = inverted_ ?
        (center_-gc)/mag(center_-gc) :
        (gc-center_)/mag(gc-center_);

    // Return gc plus twice the normal vector times the distance
    return (gc + 2.0*n*dist);
}

}

}

}
