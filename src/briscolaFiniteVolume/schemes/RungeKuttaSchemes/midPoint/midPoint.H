#ifndef midPoint_H
#define midPoint_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class midPoint
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from mesh
    midPoint(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "midPoint.C"
#endif

#endif
