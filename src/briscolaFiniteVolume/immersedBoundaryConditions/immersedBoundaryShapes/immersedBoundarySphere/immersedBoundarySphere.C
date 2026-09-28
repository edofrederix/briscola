#include "immersedBoundarySphere.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

immersedBoundarySphere::immersedBoundarySphere
(
    const dictionary& dict,
    bool inverted
)
:
    immersedBoundaryShape(dict,inverted),
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

immersedBoundarySphere::~immersedBoundarySphere()
{}

bool immersedBoundarySphere::isInside(vector point) const
{
    return (mag(center_ - point) <= radius_) != this->inverted_;
}

scalar immersedBoundarySphere::wallDistance(vector c, vector nb) const
{
    // Return -1 if the center point is not a fluid point
    // or if the neighboring point is not inside the sphere
    if (this->isInside(c) || !this->isInside(nb))
        return -1;

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

scalar immersedBoundarySphere::wallNormalDistance(vector gc) const
{
    // Return radius minus distance from center to ghost cell
    return (inverted_ ? mag(gc-center_) - radius_ : radius_ - mag(gc-center_));
}

vector immersedBoundarySphere::mirrorPoint(vector gc) const
{
    // Wall-normal distance
    scalar dist = this->wallNormalDistance(gc);

    // Wall-normal unit vector
    vector n = inverted_ ?
        (center_-gc)/max(mag(center_-gc),1e-10) :
        (gc-center_)/max(mag(gc-center_),1e-10);

    // Return gc plus twice the normal vector times the distance
    return (gc + 2.0*n*dist);
}

}

}

}
