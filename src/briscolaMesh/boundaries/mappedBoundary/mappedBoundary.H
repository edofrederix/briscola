#ifndef mappedBoundary_H
#define mappedBoundary_H

#include "patchBoundary.H"

namespace Foam
{

namespace briscola
{

class mappedBoundary
:
    public patchBoundary
{
    // Mapping offset
    vector mappingOffset_;

    // Disallow default bitwise assignment
    void operator=(const mappedBoundary&);

public:

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

    // Mapping offset
    inline const vector& mappingOffset() const
    {
        return mappingOffset_;
    }
};

}

}

#endif
