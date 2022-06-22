#include "edge.H"
#include "edgeKey.H"
#include "geometry.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(edge, 0);
defineRunTimeSelectionTable(edge, dictionary);

edge::edge
(
    const face& f,
    const label num,
    const labelBlock& vertexNums,
    const labelVector& N
)
:
    meshObject<face>(f, num),
    v_(vertexNums),
    N_(N),
    vertices_()
{
    vertices_.setSize(2);

    // Walk through vertices in a way that matches the vertexNumsInEdge array.
    // Note that v_ is squeezed.

    vertices_.set
    (
        0,
        new vertex
        (
            *this,
            vertexNumsInEdge[this->num()][0],
            vertexNums(0)
        )
    );

    vertices_.set
    (
        1,
        new vertex
        (
            *this,
            vertexNumsInEdge[this->num()][1],
            vertexNums(1)
        )
    );

    if (mag(dist()) == 0.0)
    {
        FatalErrorInFunction
            << "Zero length for " << *this
            << exit(FatalError);
    }
}

edge::edge
(
    const edge& e
)
:
    meshObject<face>(e.parentFace(), e.num()),
    v_(e.v_),
    N_(e.N_),
    vertices_(e.vertices_)
{}

edge::~edge()
{}

autoPtr<edge> edge::New
(
    const face& f,
    const label num,
    const labelBlock& vertexNums,
    const labelVector& N
)
{
    const edgeTable& edgeData = f.parentBrick().parentGeometry().edgeData();

    const edgeKey key(vertexNums(0), vertexNums(1));

    word edgeType("line");

    if (edgeData.found(key))
    {
        edgeType = word(edgeData[key].lookup("type"));
    }
    else if (edgeData.found(reverse(key)))
    {
        edgeType = word(edgeData[reverse(key)].lookup("type"));
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(edgeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown edge type " << edgeType << " for edge " << key
            << " of " << f << endl << endl
            << "Valid edge types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<edge>(cstrIter()(f, num, vertexNums, N));
}

bool edge::operator==(const edge& e) const
{
    if
    (
        (
            v0() == e.v0()
         && v1() == e.v1()
        )
     || (
            v1() == e.v0()
         && v0() == e.v1()
        )
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

dictionary edge::dict() const
{
    const edgeTable& edgeData =
        parentFace().parentBrick().parentGeometry().edgeData();

    const edgeKey key(v_(0), v_(1));

    if (edgeData.found(key))
    {
        return dictionary(edgeData[key]);
    }
    else
    {
        return dictionary();
    }
}

bool edge::operator!=(const edge& e) const
{
    return !(*this == e);
}

Ostream& operator<<(Ostream& os, const edge& e)
{
    os  << "edge " << e.num() << ": "
        << "(" << e.v()(0) << "," << e.v()(1) << ")";

    return os;
}

}

}
