#ifndef periodicBoundary_H
#define periodicBoundary_H

#include "parallelBoundary.H"

namespace Foam
{

namespace briscola
{

class periodicBoundary
:
    public parallelBoundary
{
    // Disallow default bitwise assignment
    void operator=(const periodicBoundary&);

public:

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicBoundary(const level& lvl, const dictionary& dict);

    // Copy constructor
    periodicBoundary(const periodicBoundary&);

    // Copy constructor with alternative level reference
    periodicBoundary(const periodicBoundary&, const level&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new periodicBoundary(*this));
    }

    // Clone boundary with alternative level reference
    virtual autoPtr<boundary> clone(const level& lvl) const
    {
        return autoPtr<boundary>(new periodicBoundary(*this, lvl));
    }
};

}

}

#endif
