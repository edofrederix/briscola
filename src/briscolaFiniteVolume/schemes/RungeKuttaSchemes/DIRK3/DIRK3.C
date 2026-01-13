#include "DIRK3.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

DIRK3::DIRK3(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,3)
{
    const scalar gamma(0.4358665215);
    const scalar b1(-3.0/2.0*Foam::sqr(gamma)+4.0*gamma-0.25);
    const scalar b2(3.0/2.0*Foam::sqr(gamma)-5.0*gamma+5.0/4.0);

    a_(0,0) = gamma;

    a_(1,0) = 0.5-gamma/2.0;
    a_(1,1) = gamma;

    a_(2,0) = b1;
    a_(2,1) = b2;
    a_(2,2) = gamma;

    b_ = a_;
}

}

}

}

}
