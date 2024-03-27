#include "cylinder.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

cylinder::cylinder
(
    const dictionary& dict,
    bool inverted
)
:
    shape(dict,inverted),
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

cylinder::~cylinder()
{}

bool cylinder::isInside(vector point)
{
    // Check if point is between the two planes defined by the two
    // cylinder ends, then check if the distance to the cylinder
    // axis is smaller than the radius
    if
    (
        (((point - start_) & (end_ - start_)) >= 0.0)
        && (((point - end_)   & (end_ - start_)) <= 0.0)
        && (
            (mag((point - start_) ^ (end_ - start_))
            / mag(end_ - start_)) <= radius_
        )
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

    if(!this->inverted_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

scalar cylinder::wallDistance(vector c, vector nb)
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

    // Project c onto the cylinder axis
    vector cProj = start_ + ((c-start_) & axis) * axis
        / (axis & axis);

    // Project nb onto the cylinder axis
    vector nbProj = start_ + ((nb-start_) & axis) * axis
        / (axis & axis);

    // First check if the line is aligned with the axis and
    // only intersects at the end caps of the cylinder
    if
    (
        (mag(c - cProj) == mag(nb - nbProj))
        && (mag(c-cProj) < radius_)
    )
    {
        scalar t1 = (c-start_) & C;
        scalar t2 = (c-end_) & C;

        return min(mag(t1), mag(t2));
    }

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

        // Check whether intersection points are before cylinder start
        if(((i1 - start_) & C) < 0)
        {
            t1 = ((start_-c) & C) / (D & C);
        }
        if(((i2 - start_) & C) < 0)
        {
            t2 = ((start_-c) & C) / (D & C);
        }

        // Check whether intersection points are after cylinder end
        if(((i1 - start_) & C) > mag(axis))
        {
            t1 = ((end_-c) & C) / (D & C);
        }
        if(((i2 - start_) & C) > mag(axis))
        {
            t2 = ((end_-c) & C) / (D & C);
        }

        return min(mag(t1), mag(t2));
    }
}

scalar cylinder::wallNormalDistance(vector gc)
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

vector cylinder::mirrorPoint(vector gc)
{
    scalar dist = this->wallNormalDistance(gc);

    if (dist < 0)
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
        (p-gc)/mag(p-gc) :
        (gc-p)/mag(gc-p);

    // Return gc plus twice the normal vector
    return (gc + 2.0*n*dist);
}

}

}

}
