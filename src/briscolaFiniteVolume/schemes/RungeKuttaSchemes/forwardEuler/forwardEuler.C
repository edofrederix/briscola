#include "forwardEuler.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

forwardEuler::forwardEuler(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh),
    a_(2, Zero)
{
    a_(1,0) = 1.0;
}

}

}

}

}
