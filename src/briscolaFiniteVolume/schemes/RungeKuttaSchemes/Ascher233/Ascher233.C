#include "Ascher233.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher233::Ascher233(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,4)
{
    const scalar gamma(0.5 + Foam::sqrt(3.0)/6.0);

    a_(1,0) = gamma;

    a_(2,0) = gamma - 1.0;
    a_(2,1) = 2.0*(1.0 - gamma);

    a_(3,1) = 0.5;
    a_(3,2) = 0.5;

    b_(1,1) = gamma;

    b_(2,1) = 1.0 - 2.0*gamma;
    b_(2,2) = gamma;

    b_(3,1) = 0.5;
    b_(3,2) = 0.5;
}

}

}

}

}
