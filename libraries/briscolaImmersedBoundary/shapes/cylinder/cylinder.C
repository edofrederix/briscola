#include "cylinder.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

// Constructor

cylinder::cylinder
(
    vector start,
    vector end,
    scalar radius,
    bool inverted
)
:
    start_(start),
    end_(end),
    radius_(radius),
    inverted_(inverted)
{
    if (start_ == end_)
    {
        FatalError
            << "Cylinder start and end can't coincide."
            << endl;
        FatalError.exit();
    }

    if (radius <= 0.0)
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
        if(!inverted_)
        {
            return true;
        }
        else
        {
            return false;
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

    // Normalized direction vector of the line
    vector D = (nb-c)/mag(nb-c);

    // Normalized direction vector of the cylinder axis
    vector C = (end_ - start_)/mag(end_ - start_);

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
        scalar t1 = (-y + sqrt(sqr(y)-4.0*x*z))/(2.0*x);
        scalar t2 = (-y - sqrt(sqr(y)-4.0*x*z))/(2.0*x);

        return min(mag(t1), mag(t2));
    }
}

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
