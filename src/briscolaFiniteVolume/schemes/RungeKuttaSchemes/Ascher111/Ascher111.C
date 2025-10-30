#include "Ascher111.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher111::Ascher111(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh),
    a_(2, Zero),
    b_(2, Zero)
{
    a_(1,0) = 1.0;
    b_(1,1) = 1.0;
}

}

}

}

}
