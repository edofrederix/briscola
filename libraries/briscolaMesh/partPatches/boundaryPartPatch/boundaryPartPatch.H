#ifndef boundaryPartPatch_H
#define boundaryPartPatch_H

#include "partPatch.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

class boundaryPartPatch
:
    public partPatch
{
    // Copy constructor
    boundaryPartPatch(const boundaryPartPatch&);

    // Disallow default bitwise assignment
    void operator=(const boundaryPartPatch&);

public:

    // Runtime type information
    TypeName("boundary");

    // Construct from components
    boundaryPartPatch(const mesh& msh, const dictionary& dict);
};

}

}

}

#endif
