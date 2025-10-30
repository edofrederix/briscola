#ifndef backwardEuler_H
#define backwardEuler_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class backwardEuler
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("backwardEuler");

    // Construct from mesh
    backwardEuler(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "backwardEuler.C"
#endif

#endif
