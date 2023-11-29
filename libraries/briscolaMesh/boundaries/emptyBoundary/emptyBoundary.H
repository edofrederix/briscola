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

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("empty");

    // Construct from components
    emptyBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    emptyBoundary(const emptyBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new emptyBoundary(*this));
    }

    // Boundary type number
    virtual label typeNum() const
    {
        return typeNumber;
    }
};

}

}

#endif
