#ifndef periodicPartPatch_H
#define periodicPartPatch_H

#include "parallelPartPatch.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

class periodicPartPatch
:
    public parallelPartPatch
{
    // Disallow default bitwise assignment
    void operator=(const periodicPartPatch&);

public:

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicPartPatch(const mesh& msh, const dictionary& dict);

    // Copy constructor
    periodicPartPatch(const periodicPartPatch&);

    // Clone part patch
    virtual autoPtr<partPatch> clone() const
    {
        return autoPtr<partPatch>(new periodicPartPatch(*this));
    }
};

}

}

}

#endif
