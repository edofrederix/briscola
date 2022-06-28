#ifndef parallelPartPatch_H
#define parallelPartPatch_H

#include "partPatch.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

class parallelPartPatch
:
    public partPatch
{
    // Processor number of the other side
    label neighborProcNum_;

    // Copy constructor
    parallelPartPatch(const parallelPartPatch&);

    // Disallow default bitwise assignment
    void operator=(const parallelPartPatch&);

public:

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelPartPatch(const mesh& msh, const dictionary& dict);
};

}

}

}

#endif
