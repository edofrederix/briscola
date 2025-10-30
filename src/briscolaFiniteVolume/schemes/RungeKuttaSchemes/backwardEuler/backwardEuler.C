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
    RungeKuttaScheme(fvMsh,1)
{
    a_(0,0) = 1.0;
    b_ = a_;
}

}

}

}

}
