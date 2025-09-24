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
    patchBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    patchBoundary(const patchBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new patchBoundary(*this));
    }
};

}

}

#endif
