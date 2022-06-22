#include "face.H"
#include "geometry.H"

namespace Foam
{

namespace briscola
{

void face::createEdges()
{
    edges_.setSize(4);

    label l = 0;

    // Walk through edges in a way that matches the edgeNumsInFace array. Note
    // that v_ is squeezed.

    for (label dir = 1; dir >= 0; dir--)
    {
        for (label i = 0; i < 2; i++)
        {
            labelVector N(N_);
            N[dir] = 1;

            edges_.set
            (
                l,
                edge::New
                (
                    *this,
                    edgeNumsInFace[this->num()][l],
                    v_.slice(i,dir,true),
                    N
                )
            );

            l++;
        }
    }
}

face::face
(
    const brick& b,
    const label num,
    const labelBlock& vertexNums,
    const labelVector& N
)
:
    meshObject<brick>(b, num),
    v_(vertexNums),
    N_(N),
    edges_()
{
    createEdges();

    if (area() == 0.0)
    {
        FatalErrorInFunction
            << "Zero area for " << *this
            << exit(FatalError);
    }
}

face::face
(
    const face& f
)
:
    meshObject<brick>(f.parentBrick(), f.num()),
    v_(f.v_),
    N_(f.N_),
    edges_()
{
    createEdges();
}

face::~face()
{}

bool face::operator==(const face& f) const
{
    for (label i = 0; i < 4; i++)
    {
        for (label j = 0; j < 4; j++)
        {
            if (v_(i) == f.v()(j))
            {
                break;
            }
            else if (j == 3)
            {
                return false;
            }
        }
    }

    return true;
}

bool face::operator!=(const face& f) const
{
    return !(*this == f);
}

Ostream& operator<<(Ostream& os, const face& f)
{
    os  << "face " << f.num() << ": "
        << "((" << f.v()(0) << "," << f.v()(1) << "),"
        << "(" << f.v()(2) << "," << f.v()(3) << "))";

    return os;
}

}

}
