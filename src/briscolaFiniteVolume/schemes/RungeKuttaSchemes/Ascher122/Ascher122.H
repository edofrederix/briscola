#ifndef Ascher122_H
#define Ascher122_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher122
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Ascher122");

    // Construct from mesh
    Ascher122(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher122.C"
#endif

#endif
