#ifndef DIRK2_H
#define DIRK2_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class DIRK2
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("DIRK2");

    // Construct from mesh
    DIRK2(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "DIRK2.C"
#endif

#endif
