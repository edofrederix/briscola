#ifndef mappedBoundary_H
#define mappedBoundary_H

#include "domainBoundary.H"

namespace Foam
{

namespace briscola
{

class mappedBoundary
:
    public domainBoundary
{
    // Mapping offset
    vector mappingOffset_;

    // Disallow default bitwise assignment
    void operator=(const mappedBoundary&);

public:

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("mapped");

    // Construct from components
    mappedBoundary(const mesh& msh, const dictionary& dict);

    // Copy constructor
    mappedBoundary(const mappedBoundary&);

    // Clone boundary
    virtual autoPtr<boundary> clone() const
    {
        return autoPtr<boundary>(new mappedBoundary(*this));
    }

    // Boundary type number
    virtual label typeNum() const
    {
        return typeNumber;
    }

    // Mapping offset
    inline const vector& mappingOffset() const
    {
        return mappingOffset_;
    }
};

}

}

#endif
