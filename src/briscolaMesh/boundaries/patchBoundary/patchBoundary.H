#ifndef patchBoundary_H
#define patchBoundary_H

#include "boundary.H"

namespace Foam
{

namespace briscola
{

class patchBoundary
:
    public boundary
{
    // Disallow default bitwise assignment
    void operator=(const patchBoundary&);

public:

    // Runtime type information
    TypeName("patch");

    // Construct from components
    patchBoundary(const level& lvl, const dictionary& dict);

    // Copy constructor
    patchBoundary(const patchBoundary&);

    // Copy constructor with alternative level reference
    patchBoundary(const patchBoundary&, const level&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new patchBoundary(*this));
    }

    // Clone boundary with alternative level reference
    virtual autoPtr<boundary> clone(const level& lvl) const
    {
        return autoPtr<boundary>(new patchBoundary(*this, lvl));
    }
};

}

}

#endif
