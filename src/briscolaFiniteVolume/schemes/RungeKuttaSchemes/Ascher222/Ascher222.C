#include "Ascher222.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher222::Ascher222(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh),
    a_(3, Zero),
    b_(3, Zero)
{
    const scalar gamma(1.0 - Foam::sqrt(2.0)/2.0);
    const scalar delta(1.0 - 0.5/gamma);

    a_(1,0) = gamma;

    a_(2,0) = delta;
    a_(2,1) = 1.0 - delta;

    b_(1,1) = gamma;

    b_(2,1) = 1.0 - gamma;
    b_(2,2) = gamma;
}

}

}

}

}
