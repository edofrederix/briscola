#include "RK4CN.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

RK4CN::RK4CN(const fvMesh& fvMsh)
:
    RK4(fvMsh)
{
    b_.setSize(5, scalarList(5, 0.0));

    b_[4][0] = 1.0/2.0;
    b_[4][4] = 1.0/2.0;
}

}

}

}

}
