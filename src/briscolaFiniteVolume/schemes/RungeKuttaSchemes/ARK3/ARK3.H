#ifndef ARK3_H
#define ARK3_H

#include "RungeKuttaScheme.H"

/**

@brief Accelerated Runge-Kutta 3 (ARK3) scheme

See Udwadia, F. E., & Farahani, A. (2008). Accelerated Runge‐Kutta Methods.
Discrete Dynamics in Nature and Society, 2008(1), 790619.

*/

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class ARK3
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("ARK3");

    // Construct from mesh
    ARK3(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "ARK3.C"
#endif

#endif
