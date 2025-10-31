#include "CNAB.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

CNAB::CNAB(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh,2,2)
{
    a_(1,0) = 1.5;
    a_(1,2) = -0.5;
    b_(1,0) = 0.5;
    b_(1,1) = 0.5;
}

}

}

}

}
