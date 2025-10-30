#ifndef Ascher121_H
#define Ascher121_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher121
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("Ascher121");

    // Construct from mesh
    Ascher121(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher121.C"
#endif

#endif
