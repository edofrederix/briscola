#ifndef Ascher111_H
#define Ascher111_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher111
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Ascher111");

    // Construct from mesh
    Ascher111(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher111.C"
#endif

#endif
