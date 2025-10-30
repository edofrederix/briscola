#ifndef RK4_H
#define RK4_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class RK4
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("RK4");

    // Construct from mesh
    RK4(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "RK4.C"
#endif

#endif
