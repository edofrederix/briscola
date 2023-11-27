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

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    periodicBoundary(const periodicBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new periodicBoundary(*this));
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
