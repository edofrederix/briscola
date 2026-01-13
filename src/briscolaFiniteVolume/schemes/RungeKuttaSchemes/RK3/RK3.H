#ifndef RK3_H
#define RK3_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class RK3
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("RK3");

    // Construct from mesh
    RK3(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "RK3.C"
#endif

#endif
