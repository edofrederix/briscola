#ifndef emptyBoundary_H
#define emptyBoundary_H

#include "boundary.H"

namespace Foam
{

namespace briscola
{

class emptyBoundary
:
    public boundary
{
    // Disallow default bitwise assignment
    void operator=(const emptyBoundary&);

public:

    // Runtime type information
    TypeName("empty");

    // Construct from components
    emptyBoundary(const level& lvl, const dictionary& dict);

    // Copy constructor
    emptyBoundary(const emptyBoundary&);

    // Copy constructor with alternative level reference
    emptyBoundary(const emptyBoundary&, const level&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new emptyBoundary(*this));
    }

    // Clone boundary with alternative level reference
    virtual autoPtr<boundary> clone(const level& lvl) const
    {
        return autoPtr<boundary>(new emptyBoundary(*this, lvl));
    }
};

}

}

#endif
