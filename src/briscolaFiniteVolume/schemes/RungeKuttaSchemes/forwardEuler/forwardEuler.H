#ifndef forwardEuler_H
#define forwardEuler_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class forwardEuler
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("forwardEuler");

    // Construct from mesh
    forwardEuler(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "forwardEuler.C"
#endif

#endif
