#include "RK3CN.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

RK3CN::RK3CN(const fvMesh& fvMsh)
:
    RK3(fvMsh)
{
    b_.setSize(4, scalarList(4, 0.0));

    b_[3][0] = 1.0/2.0;
    b_[3][3] = 1.0/2.0;
}

}

}

}

}
