#include "Ascher111.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher111::Ascher111(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(2, scalarList(2, 0.0));
    b_.setSize(2, scalarList(2, 0.0));

    a_[1][0] = 1.0;
    b_[1][1] = 1.0;
}

}

}

}

}
