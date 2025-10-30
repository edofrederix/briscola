#ifndef DIRK4_H
#define DIRK4_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class DIRK4
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("DIRK4");

    // Construct from mesh
    DIRK4(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "DIRK4.C"
#endif

#endif
