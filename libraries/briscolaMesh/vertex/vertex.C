#include "vertex.H"
#include "geometry.H"

namespace Foam
{

namespace briscola
{

vertex::vertex(const edge& e, const label num, const label j)
:
    meshObject<edge>(e, num),
    vector
    (
        e
       .parentFace()
       .parentBrick()
       .parentGeometry()
       .vertexData()[j]
    ),
    j_(j)
{}

vertex::vertex(const vertex& v)
:
    meshObject<edge>(v.parentEdge(), v.num()),
    vector(v),
    j_(v.j_)
{}

vertex::vertex(const vertex& v, const edge& e)
:
    meshObject<edge>(e, v.num()),
    vector(v),
    j_(v.j_)
{}

vertex::~vertex()
{}

void vertex::operator=(const vertex& v)
{
    this->v() = v.v();
}

Ostream& operator<<(Ostream& os, const vertex& v)
{
    os  << "vertex " << v.num() << ": " << v.j();

    return os;
}

}

}
