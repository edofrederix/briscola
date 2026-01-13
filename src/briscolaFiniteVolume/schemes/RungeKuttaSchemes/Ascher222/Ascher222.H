#ifndef Ascher222_H
#define Ascher222_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher222
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Ascher222");

    // Construct from mesh
    Ascher222(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher222.C"
#endif

#endif
