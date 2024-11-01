#include "forwardEuler.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

forwardEuler::forwardEuler(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(2, scalarList(2, 0.0));

    a_[1][0] = 1.0;
}

}

}

}

}
