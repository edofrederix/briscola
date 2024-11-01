#ifndef CrankNicolson_H
#define CrankNicolson_H

#include "theta.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class CrankNicolson
:
    public theta
{
public:

    // Runtime type information
    TypeName("CrankNicolson");

    // Construct from mesh
    CrankNicolson(const fvMesh& fvMsh);

    // Destructor
    ~CrankNicolson()
    {}
};

}

}

}

}

#ifdef NoRepository
    #include "CrankNicolson.C"
#endif

#endif
