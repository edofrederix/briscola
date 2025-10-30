#ifndef DIRK3_H
#define DIRK3_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class DIRK3
:
    public RungeKuttaScheme
{
public:

    // Runtime type information
    TypeName("DIRK3");

    // Construct from mesh
    DIRK3(const fvMesh& fvMsh);
};

}

}

}

}

#ifdef NoRepository
    #include "DIRK3.C"
#endif

#endif
