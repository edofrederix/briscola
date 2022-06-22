#ifndef brickLinks_H
#define brickLinks_H

#include "brick.H"
#include "brickFaceLink.H"
#include "brickEdgeLink.H"
#include "brickVertexLink.H"

namespace Foam
{

namespace briscola
{

class brickLinks
{
    // Reference to the brick
    const brick& b_;

    // Face links
    PtrList<brickFaceLink> faceLinks_;

    // Edge links
    PtrList<brickEdgeLink> edgeLinks_;

    // Vertex links
    PtrList<brickVertexLink> vertexLinks_;

public:

    // Construct from bricks
    brickLinks(const brick&);

    // Destructor
    virtual ~brickLinks();

    // Access to the face links
    PtrList<brickFaceLink>& faceLinks() { return faceLinks_; }
    const PtrList<brickFaceLink>& faceLinks() const { return faceLinks_; }

    // Access to the edge links
    PtrList<brickEdgeLink>& edgeLinks() { return edgeLinks_; }
    const PtrList<brickEdgeLink>& edgeLinks() const { return edgeLinks_; }

    // Access to the vertex links
    PtrList<brickVertexLink>& vertexLinks() { return vertexLinks_; }
    const PtrList<brickVertexLink>& vertexLinks() const { return vertexLinks_; }

    // Find face number of link with given brick. Return -1 if no such link
    // exists.
    label faceLinkNum(const label) const;

    // Return reference to face link with given brick number
    const brickFaceLink& getFaceLink(const label) const;

    // Find edge number of link with given brick. Return -1 if no such link
    // exists.
    label edgeLinkNum(const label) const;

    // Return reference to edge link with given brick number
    const brickEdgeLink& getEdgeLink(const label) const;

    // Find vertex number of link with given brick. Return -1 if no such link
    // exists.
    label vertexLinkNum(const label) const;

    // Return reference to vertex link with given brick number
    const brickVertexLink& getVertexLink(const label) const;
};

}

}

#endif
