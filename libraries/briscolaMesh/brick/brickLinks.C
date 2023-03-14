#include "brickLinks.H"

namespace Foam
{

namespace briscola
{

brickLinks::brickLinks(const brick& b)
:
    b_(b),
    faceLinks_(6),
    edgeLinks_(12),
    vertexLinks_(8)
{}

brickLinks::brickLinks(const brickLinks& links)
:
    b_(links.b_),
    faceLinks_(links.faceLinks_),
    edgeLinks_(links.edgeLinks_),
    vertexLinks_(links.vertexLinks_)
{}

brickLinks::~brickLinks()
{}

label brickLinks::faceLinkNum(const label bricki) const
{
    forAll(faceLinks_, facei)
    {
        if
        (
            faceLinks_.set(facei)
         && faceLinks_[facei].f1().parentBrick().num() == bricki
        )
        {
            return facei;
        }
    }

    return -1;
}

const brickFaceLink& brickLinks::getFaceLink(const label bricki) const
{
    const label facei = faceLinkNum(bricki);

    if (facei == -1)
    {
        FatalErrorInFunction
            << "Face link requested for brick " << b_.num() << " to "
            << bricki << " but it does not exist." << endl
            << abort(FatalError);
    }

    return faceLinks_[facei];
}

label brickLinks::edgeLinkNum(const label bricki) const
{
    forAll(edgeLinks_, edgei)
    {
        if
        (
            edgeLinks_.set(edgei)
         && edgeLinks_[edgei].e1().parentFace().parentBrick().num()
         == bricki
        )
        {
            return edgei;
        }
    }

    return -1;
}

const brickEdgeLink& brickLinks::getEdgeLink(const label bricki) const
{
    const label edgei = edgeLinkNum(bricki);

    if (edgei == -1)
    {
        FatalErrorInFunction
            << "Edge link requested for brick " << b_.num() << " to "
            << bricki << " but it does not exist." << endl
            << abort(FatalError);
    }

    return edgeLinks_[edgei];
}

label brickLinks::vertexLinkNum(const label bricki) const
{
    forAll(vertexLinks_, vertexi)
    {
        if
        (
            vertexLinks_.set(vertexi)
         && vertexLinks_[vertexi].v1().parentEdge().parentFace().parentBrick().num()
         == bricki
        )
        {
            return vertexi;
        }
    }

    return -1;
}

const brickVertexLink& brickLinks::getVertexLink(const label bricki) const
{
    const label vertexi = vertexLinkNum(bricki);

    if (vertexi == -1)
    {
        FatalErrorInFunction
            << "Vertex link requested for brick " << b_.num() << " to "
            << bricki << " but it does not exist." << endl
            << abort(FatalError);
    }

    return vertexLinks_[vertexi];
}

}

}
