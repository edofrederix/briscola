#include "brickLink.H"

namespace Foam
{

namespace briscola
{

brickLink::brickLink
(
    const brick& b0,
    const brick& b1,
    const labelVector offset,
    const labelTensor T,
    const bool periodic
)
:
    b0_(b0),
    b1_(b1),
    offset_(offset),
    T_(T),
    periodic_(periodic)
{}

brickLink::brickLink(const brickLink& link)
:
    b0_(link.b0_),
    b1_(link.b1_),
    offset_(link.offset_),
    T_(link.T_),
    periodic_(link.periodic_)
{}

brickLink::brickLink(const brickFaceLink& link)
:
    b0_(link.f0().parentBrick()),
    b1_(link.f1().parentBrick()),
    offset_(faceOffsets[link.f0().num()]),
    T_(link.T()),
    periodic_(link.periodic())
{}

brickLink::brickLink(const brickEdgeLink& link)
:
    b0_(link.e0().parentFace().parentBrick()),
    b1_(link.e1().parentFace().parentBrick()),
    offset_(edgeOffsets[link.e0().num()]),
    T_(link.T()),
    periodic_(link.periodic())
{}

brickLink::brickLink(const brickVertexLink& link)
:
    b0_(link.v0().parentEdge().parentFace().parentBrick()),
    b1_(link.v1().parentEdge().parentFace().parentBrick()),
    offset_(vertexOffsets[link.v0().num()]),
    T_(link.T()),
    periodic_(link.periodic())
{}

brickLink::~brickLink()
{}

}

}
