#include "brickVertexLink.H"
#include "brickTopology.H"

namespace Foam
{

namespace briscola
{

brickVertexLink::brickVertexLink
(
    const vertex& v0,
    const vertex& v1,
    const labelTensor& T,
    const bool periodic
)
:
    v0_(v0),
    v1_(v1),
    T_(T),
    periodic_(periodic)
{}

brickVertexLink::~brickVertexLink()
{}

}

}
