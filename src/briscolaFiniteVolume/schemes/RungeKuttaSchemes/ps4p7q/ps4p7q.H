#ifndef ps4p7q_H
#define ps4p7q_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class ps4p7q
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("ps4p7q");

    // Construct from mesh
    ps4p7q(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "ps4p7q.C"
#endif

#endif
