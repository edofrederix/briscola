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

    // Master bool. This part patch is master if its face, edge or vertex number
    // is even and the neighboring one is odd. If they are both odd or even, the
    // lower processor number is master.
    bool master_;

    // Copy constructor
    parallelPartPatch(const parallelPartPatch&);

    // Disallow default bitwise assignment
    void operator=(const parallelPartPatch&);

public:

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelPartPatch(const mesh& msh, const dictionary& dict);

    // Master/slave bools
    virtual bool master() const { return master_; }
    virtual bool slave() const { return !master_; }
};

}

}

}

#endif
