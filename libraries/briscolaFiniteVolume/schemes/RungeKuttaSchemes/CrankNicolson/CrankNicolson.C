#include "CrankNicolson.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

CrankNicolson::CrankNicolson(const fvMesh& fvMsh)
:
    theta(fvMsh, 0.5)
{}

}

}

}

}
