#ifndef arcEdge_H
#define arcEdge_H

#include "edge.H"
#include "constants.H"

namespace Foam
{

namespace briscola
{

namespace edges
{

class arcEdge
:
    public edge
{
    // Center of the arc
    vector center_;

    // Angle of the arc
    scalar angle_;

    // Unit axis of rotation
    vector axis_;

    // Radius of the arc
    scalar radius_;

    // Unit start vector relative to center
    vector start_;

    // Copy constructor
    arcEdge(const arcEdge&);

    // Disallow default bitwise assignment
    void operator=(const arcEdge&);

public:

    // Runtime type information
    TypeName("arc");

    // Construct from components
    arcEdge(const face&, const label, const labelBlock&, const labelVector&);

    // Length of the edge
    virtual scalar length() const
    {
        return angle_*mag(start_);
    }

    // Parametrized coordinate
    virtual tmp<vectorField> operator()(const scalarField&) const;
};

}

}

}

#endif
