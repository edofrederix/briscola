#ifndef emptyPartPatch_H
#define emptyPartPatch_H

#include "partPatch.H"

namespace Foam
{

namespace briscola
{

class emptyPartPatch
:
    public partPatch
{
    // Disallow default bitwise assignment
    void operator=(const emptyPartPatch&);

public:

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("empty");

    // Construct from components
    emptyPartPatch(const mesh& msh, const dictionary& dict);

    // Copy constructor
    emptyPartPatch(const emptyPartPatch&);

    // Clone part patch
    virtual autoPtr<partPatch> clone() const
    {
        return autoPtr<partPatch>(new emptyPartPatch(*this));
    }

    // Part patch type number
    virtual label typeNum() const
    {
        return typeNumber;
    }
};

}

}

#endif
