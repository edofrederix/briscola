#include "RK3.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

RK3::RK3(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,4)
{
    a_(1,0) = 0.5;

    a_(2,0) = -1.0;
    a_(2,1) = 2.0;

    a_(3,0) = 1.0/6.0;
    a_(3,1) = 2.0/3.0;
    a_(3,2) = 1.0/6.0;

    b_ = a_;
}

}

}

}

}
