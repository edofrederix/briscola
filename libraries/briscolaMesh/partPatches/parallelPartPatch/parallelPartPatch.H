#ifndef parallelPartPatch_H
#define parallelPartPatch_H

#include "partPatch.H"

namespace Foam
{

namespace briscola
{

class parallelPartPatch
:
    public partPatch
{
    // Processor number of the other side
    label neighborProcNum_;

    // Disallow default bitwise assignment
    void operator=(const parallelPartPatch&);

public:

    // Type number
    const static label typeNumber;

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelPartPatch(const mesh& msh, const dictionary& dict);

    // Copy constructor
    parallelPartPatch(const parallelPartPatch&);

    // Clone part patch
    virtual autoPtr<partPatch> clone() const
    {
        return autoPtr<partPatch>(new parallelPartPatch(*this));
    }

    // Part patch type number
    virtual label typeNum() const
    {
        return typeNumber;
    }

    // Neighbor processor number
    virtual label neighborProcNum() const { return neighborProcNum_; }
};

}

}

#endif
