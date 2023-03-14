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

    // Copy constructors with new face
    lineEdge(const lineEdge&, const face&);

    // Copy constructor from edge, because lineEdge contains no data
    lineEdge(const edge&);

    // Clone edge
    virtual autoPtr<edge> clone() const
    {
        return autoPtr<edge>(new lineEdge(*this));
    }

    // Clone edge with new face
    virtual autoPtr<edge> clone(const face& f) const
    {
        return autoPtr<edge>(new lineEdge(*this, f));
    }

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
