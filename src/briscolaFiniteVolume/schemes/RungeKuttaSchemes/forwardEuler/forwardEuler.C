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
    RungeKuttaScheme(fvMsh,2)
{
    a_(1,0) = 1.0;
    b_ = a_;
}

}

}

}

}
