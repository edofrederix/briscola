#include "ARK3.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

ARK3::ARK3(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,3,2)
{
    a_(1,0) = 5.0/12.0;

    a_(2,0) = 0.5;
    a_(2,1) = 1.0;
    a_(2,3) = 0.5;
    a_(2,4) = -1.0;

    b_ = a_;
}

}

}

}

}
