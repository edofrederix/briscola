#include "Ascher232.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

Ascher232::Ascher232(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(4, scalarList(4, 0.0));
    b_.setSize(4, scalarList(4, 0.0));

    const scalar gamma(1.0-Foam::sqrt(2.0)/2.0);
    const scalar delta(-2.0*Foam::sqrt(2.0)/3.0);

    a_[1][0] = gamma;

    a_[2][0] = delta;
    a_[2][1] = 1.0 - delta;

    a_[3][1] = 1.0 - gamma;
    a_[3][2] = gamma;

    b_[1][1] = gamma;

    b_[2][1] = 1.0 - gamma;
    b_[2][2] = gamma;

    b_[3][1] = 1.0 - gamma;
    b_[3][2] = gamma;
}

}

}

}

}
