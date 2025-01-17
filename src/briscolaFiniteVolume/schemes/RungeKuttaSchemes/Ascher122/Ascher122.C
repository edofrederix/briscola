#include "Ascher122.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher122::Ascher122(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(3, scalarList(3, 0.0));
    b_.setSize(3, scalarList(3, 0.0));

    a_[1][0] = 0.5;
    a_[2][1] = 1.0;

    b_[1][1] = 0.5;
    b_[2][1] = 1.0;
}

}

}

}

}
