#ifndef AB2_H
#define AB2_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class AB2
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("AB2");

    // Construct from mesh
    AB2(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "AB2.C"
#endif

#endif
