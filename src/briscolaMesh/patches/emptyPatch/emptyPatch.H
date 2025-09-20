#ifndef emptyPatch_H
#define emptyPatch_H

#include "patch.H"

namespace Foam
{

namespace briscola
{

class emptyPatch
:
    public patch
{
public:

    // Runtime type information
    TypeName("empty");

    // Construct from dictionary
    emptyPatch(const geometry&, const label, const word, const dictionary&);

    // Copy constructor
    emptyPatch(const emptyPatch&);

    // Clone this empty patch
    inline autoPtr<emptyPatch> clone() const
    {
        return autoPtr<emptyPatch>(new emptyPatch(*this));
    }

    // Destructor
    virtual ~emptyPatch();
};

}

}

#endif
