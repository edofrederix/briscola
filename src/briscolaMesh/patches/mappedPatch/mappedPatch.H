#ifndef mappedPatch_H
#define mappedPatch_H

#include "patch.H"

namespace Foam
{

namespace briscola
{

class mappedPatch
:
    public patch
{
public:

    // Runtime type information
    TypeName("mapped");

    // Construct from dictionary
    mappedPatch(const geometry&, const label, const word, const dictionary&);

    // Copy constructor
    mappedPatch(const mappedPatch&);

    // Clone this mapped patch
    inline autoPtr<mappedPatch> clone() const
    {
        return autoPtr<mappedPatch>(new mappedPatch(*this));
    }

    // Destructor
    virtual ~mappedPatch();
};

}

}

#endif
