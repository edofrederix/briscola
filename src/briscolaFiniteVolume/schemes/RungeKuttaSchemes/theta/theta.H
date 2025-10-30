#ifndef theta_H
#define theta_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class theta
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("theta");

    // Construct from mesh
    theta(const fvMesh& fvMsh);

    // Construct from mesh and theta value
    theta
    (
        const fvMesh& fvMsh,
        const scalar theta
    );
};

}

}

}

}

#ifdef NoRepository
    #include "theta.C"
#endif

#endif
