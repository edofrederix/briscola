#include "theta.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

theta::theta(const fvMesh& fvMsh)
:
    theta(fvMsh, readScalar(this->dict_.lookup("theta")))
{}

theta::theta(const fvMesh& fvMsh, const scalar t)
:
    RungeKuttaScheme(fvMsh, 2)
{
    a_(1,0) = (1.0 - t);
    a_(1,1) = t;
    b_ = a_;
}

}

}

}

}
