#ifndef lineEdge_H
#define lineEdge_H

#include "edge.H"

namespace Foam
{

namespace briscola
{

class lineEdge
:
    public edge
{
    // Disallow default bitwise assignment
    void operator=(const lineEdge&);

public:

    // Runtime type information
    TypeName("line");

    // Construct from components
    lineEdge(const face&, const label, const labelBlock&, const labelVector&);

    // Copy constructors
    lineEdge(const lineEdge&);

    // Copy constructor from edge, because lineEdge contains no data
    lineEdge(const edge&);

    // Length of the edge
    virtual scalar length() const
    {
        return mag(dist());
    }

    // Parametrized coordinates
    virtual tmp<vectorField> operator()(const scalarField& f) const
    {
        return v0()*(1.0-f) + v1()*f;
    }
};

}

}

#endif
