#include "Ascher121.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher121::Ascher121(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(3, scalarList(3, 0.0));
    b_.setSize(3, scalarList(3, 0.0));

    a_[1][0] = 1.0;
    a_[2][1] = 1.0;

    b_[1][1] = 1.0;
    b_[2][1] = 1.0;
}

}

}

}

}
