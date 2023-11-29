#ifndef domainBoundary_H
#define domainBoundary_H

#include "boundary.H"

namespace Foam
{

namespace briscola
{

class domainBoundary
:
    public boundary
{
    // Disallow default bitwise assignment
    void operator=(const domainBoundary&);

public:

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("boundary");

    // Construct from components
    domainBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    domainBoundary(const domainBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new domainBoundary(*this));
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
