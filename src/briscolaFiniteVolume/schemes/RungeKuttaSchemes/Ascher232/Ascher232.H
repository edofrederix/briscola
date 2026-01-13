#ifndef Ascher232_H
#define Ascher232_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher232
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Ascher232");

    // Construct from mesh
    Ascher232(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher232.C"
#endif

#endif
