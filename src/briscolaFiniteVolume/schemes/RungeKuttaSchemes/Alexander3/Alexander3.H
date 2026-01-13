#ifndef Alexander3_H
#define Alexander3_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Alexander3
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Alexander3");

    // Construct from mesh
    Alexander3(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Alexander3.C"
#endif

#endif
