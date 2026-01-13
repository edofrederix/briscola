#include "CNABM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

CNABM::CNABM(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,2,2)
{
    a_(1,0) = 1.5;
    a_(1,2) = -0.5;

    b_(1,0) = 6.0/16.0;
    b_(1,1) = 9.0/16.0;
    b_(1,2) = 1.0/16.0;
}

}

}

}

}
