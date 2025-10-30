#ifndef ps3p5q_H
#define ps3p5q_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class ps3p5q
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("ps3p5q");

    // Construct from mesh
    ps3p5q(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "ps3p5q.C"
#endif

#endif
