#include "DIRK2.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

DIRK2::DIRK2(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,2)
{
    const scalar gamma(1.0-Foam::sqrt(2.0)/2.0);

    a_(0,0) = gamma;

    a_(1,0) = 1.0-gamma;
    a_(1,1) = gamma;

    b_ = a_;
}

}

}

}

}
