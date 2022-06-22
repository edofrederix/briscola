#ifndef periodicPartPatch_H
#define periodicPartPatch_H

#include "partPatch.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

class periodicPartPatch
:
    public partPatch
{
    // Copy constructor
    periodicPartPatch(const periodicPartPatch&);

    // Disallow default bitwise assignment
    void operator=(const periodicPartPatch&);

public:

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicPartPatch(const mesh& msh, const dictionary& dict);
};

}

}

}

#endif
