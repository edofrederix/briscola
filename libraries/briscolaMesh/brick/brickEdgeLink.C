#include "brickEdgeLink.H"
#include "brickTopology.H"

namespace Foam
{

namespace briscola
{

brickEdgeLink::brickEdgeLink
(
    const edge& e0,
    const edge& e1,
    const labelTensor& T,
    const bool periodic
)
:
    e0_(e0),
    e1_(e1),
    T_(T),
    periodic_(periodic)
{}

brickEdgeLink::~brickEdgeLink()
{}

}

}
