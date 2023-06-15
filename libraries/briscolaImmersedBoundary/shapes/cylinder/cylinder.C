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

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
