#ifndef Hairer5_H
#define Hairer5_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Hairer5
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Hairer5");

    // Construct from mesh
    Hairer5(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Hairer5.C"
#endif

#endif
