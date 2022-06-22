#ifndef manualDecomp_H
#define manualDecomp_H

#include "decomposition.H"

namespace Foam
{

namespace briscola
{

namespace decompositions
{

class manualDecomp
:
    public decomposition
{
    // Copy constructor
    manualDecomp(const manualDecomp&);

    // Disallow default bitwise assignment
    void operator=(const manualDecomp&);

public:

    // Runtime type information
    TypeName("manual");

    // Construct from mesh
    manualDecomp(mesh&);

};

}

}

}

#endif
