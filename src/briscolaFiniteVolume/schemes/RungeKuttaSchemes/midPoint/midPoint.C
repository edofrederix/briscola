#include "midPoint.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

midPoint::midPoint(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh),
    a_(2, Zero)
{
    a_(0,0) = 0.5;
    a_(1,0) = 1.0;
}

}

}

}

}
