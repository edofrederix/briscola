#include "backwardEuler.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

backwardEuler::backwardEuler(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh),
    a_(1, Zero)
{
    a_(0,0) = 1.0;
}

}

}

}

}
