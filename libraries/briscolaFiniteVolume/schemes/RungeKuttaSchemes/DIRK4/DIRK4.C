#include "DIRK4.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

DIRK4::DIRK4(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(4, scalarList(4, 0.0));

    a_[0][0] = 0.5;

    a_[1][0] = 1.0/6.0;
    a_[1][1] = 0.5;

    a_[2][0] = -0.5;
    a_[2][1] = 0.5;
    a_[2][2] = 0.5;

    a_[3][0] = 3.0/2.0;
    a_[3][1] = -3.0/2.0;
    a_[3][2] = 0.5;
    a_[3][3] = 0.5;
}

}

}

}

}
