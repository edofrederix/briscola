#include "ps3p5q.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

ps3p5q::ps3p5q(const fvMesh& fvMsh)
:
    RungeKuttaScheme(fvMsh)
{
    a_.setSize(5, scalarList(5, 0.0));

    const scalar theta(readScalar(this->dict_.lookup("theta")));

    a_[1][0] = (theta-1)/(4*theta-3);

    const scalar a3 = (2*theta-1)*(4*theta-3)/(2*(theta-1));
    a_[2][0] = theta-a3;
    a_[2][1] = a3;

    a_[3][0] = -(2*theta-1)*(2*theta-1)/(2*(theta-1)*(4*theta-3));
    a_[3][1] = (6*theta*theta-8*theta+3)/(2*(theta-1)*(2*theta-1));
    a_[3][2] = (theta-1)/(2*theta-1)/(4*theta-3);

    a_[4][0] = 1/(12*(theta-1));
    a_[4][1] = (4*theta-3)*(4*theta-3)/(12*(theta-1)*(2*theta-1));
    a_[4][2] = -a_[4][0]/(2*theta-1);
    a_[4][3] = (4*theta-3)*a_[4][0];
}

}

}

}

}
