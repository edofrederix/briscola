#include "RK4.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

RK4::RK4(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,5)
{
    a_(1,0) = 0.5;

    a_(2,1) = 0.5;

    a_(3,2) = 1.0;

    a_(4,0) = 1.0/6.0;
    a_(4,1) = 1.0/3.0;
    a_(4,2) = 1.0/3.0;
    a_(4,3) = 1.0/6.0;

    b_ = a_;
}

}

}

}

}
