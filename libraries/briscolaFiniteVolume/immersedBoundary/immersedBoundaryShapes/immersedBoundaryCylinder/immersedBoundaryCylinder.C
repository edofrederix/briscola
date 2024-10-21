#include "immersedBoundaryCylinder.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

immersedBoundaryCylinder::immersedBoundaryCylinder
(
    const dictionary& dict,
    bool inverted
)
:
    immersedBoundaryShape(dict,inverted),
    start_(vector(dict.lookup("start"))),
    end_(vector(dict.lookup("end"))),
    radius_(readScalar(dict.lookup("radius")))
{
    if (start_ == end_)
    {
        FatalError
            << "Cylinder start and end can't coincide."
            << endl;
        FatalError.exit();
    }

    if (radius_ <= 0.0)
    {
        FatalError
            << "Cylinder radius has to be a positive value."
            << endl;
        FatalError.exit();
    }
}

// Destructor

immersedBoundaryCylinder::~immersedBoundaryCylinder()
{}

bool immersedBoundaryCylinder::isInside(vector point) const
{
    // Check if the distance to the cylinder
    // axis is smaller than the radius
    if
    (
        (mag((point - start_) ^ (end_ - start_))
        / mag(end_ - start_)) <= radius_
    )
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

scalar immersedBoundaryCylinder::wallDistance(vector c, vector nb) const
{
    // Return -1 if the center point is not a fluid point
    // or if the neighboring point is not inside the cylinder
    if (this->isInside(c))
    {
        return -1;
    }

    if (!this->isInside(nb))
    {
        return -1;
    }

    // Cylinder axis
    vector axis = end_ - start_;

    // Normalized direction vector of the line
    vector D = (nb-c)/mag(nb-c);

    // Normalized direction vector of the cylinder axis
    vector C = axis/mag(axis);

    // Vector from line origin to cylinder origin
    vector w = (c - start_);

    // We define the line through c and nb as L(t)=c+t*D
    // and solve the quadratic system xt^2 + yt + z = 0
    scalar x = (D & D) - sqr(D & C);
    scalar y = 2.0 * ((D & w) - (D & C)*(w & C));
    scalar z = (w & w) - sqr(w & C) - sqr(radius_);

    if ((sqr(y)-4.0*x*z) < 0)
    {
        return -1;
    }
    else if
    (
           ((sqr(y)-4.0*x*z) == 0)
        && (mag(D & C) == 1) // Line is on cylinder surface
    )
    {
        return -1;
    }
    else
    {
        // Intersection distances with infinite cylinder
        scalar t1 = (-y + sqrt(sqr(y)-4.0*x*z))/(2.0*x);
        scalar t2 = (-y - sqrt(sqr(y)-4.0*x*z))/(2.0*x);

        // Intersection points
        vector i1 = c + t1 * D;
        vector i2 = c + t2 * D;

        // Select intersection point between c and nb
        if (mag(c-i1)+mag(nb-i1) < mag(c-i2)+mag(nb-i2))
        {
            return mag(t1);
        }
        else
        {
            return mag(t2);
        }
    }
}

scalar immersedBoundaryCylinder::wallNormalDistance(vector gc) const
{
    // Return -1 if the point is outside of the cylinder
    if (!this->isInside(gc))
    {
        return -1;
    }

    // Normalized cylinder axis vector
    vector axis = (end_-start_)/mag(end_-start_);

    // GC vector
    vector g = (gc-start_);

    // Return radius minus distance from axis to ghost cell
    return (inverted_ ? mag(axis ^ g) - radius_ : radius_ - mag(axis ^ g));
}

vector immersedBoundaryCylinder::mirrorPoint(vector gc) const
{
    scalar dist = this->wallNormalDistance(gc);

    if (dist <= 0 || dist == radius_)
    {
        // This could be a problem if gc is exactly on the IB
        return gc;
    }

    // Normalized cylinder axis vector
    vector axis = (end_-start_)/mag(end_-start_);

    // GC vector
    vector g = (gc-start_);

    // Intersection of axis and wall-normal line through gc
    vector p = start_ + axis * (axis & g);

    // Wall-normal unit vector
    vector n = inverted_ ?
        (p-gc)/max(mag(p-gc),1e-10) :
        (gc-p)/max(mag(gc-p),1e-10);

    // Return gc plus twice the normal vector
    return (gc + 2.0*n*dist);
}

}

}

}
