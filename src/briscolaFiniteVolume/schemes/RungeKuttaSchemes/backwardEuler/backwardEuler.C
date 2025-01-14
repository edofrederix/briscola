#include "backwardEuler.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

backwardEuler::backwardEuler(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(1, scalarList(1, 0.0));

    a_[0][0] = 1.0;
}

}

}

}

}
