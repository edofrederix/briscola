#ifndef brickVertexLink_H
#define brickVertexLink_H

#include "vertex.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "labelList.H"
#include "utils.H"

namespace Foam
{

namespace briscola
{

class brickTopology;
class brickLink;

class brickVertexLink
{
    // Reference to the parent vertex
    const vertex& v0_;

    // Reference to the daughter vertex
    const vertex& v1_;

    // Transformation tensor of the daughter brick to the parent brick
    labelTensor T_;

    // Periodic link
    bool periodic_;

public:

    // Construct from parent, daughter and transformation tensor. Optionally set
    // as periodic link.
    brickVertexLink
    (
        const vertex&,
        const vertex&,
        const labelTensor&,
        const bool = false
    );

    // Copy constructor
    inline brickVertexLink(const brickVertexLink& link)
    :
        v0_(link.v0_),
        v1_(link.v1_),
        T_(link.T_),
        periodic_(link.periodic_)
    {}

    // Clone link
    inline autoPtr<brickVertexLink> clone() const
    {
        return autoPtr<brickVertexLink>(new brickVertexLink(*this));
    }

    // Destructor
    virtual ~brickVertexLink();

    // Access to the parent vertex
    inline const vertex& v0() const { return v0_; }

    // Access to the daughter vertex
    inline const vertex& v1() const { return v1_; }

    // Calculate the offset
    inline labelVector offset() const
    {
        return vertexOffsets[v0_.num()];
    }

    // Access to the transformation tensor
    inline const labelTensor& T() const { return T_; }

    // Periodic link
    inline bool periodic() const { return periodic_; }

    // Generate brick link
    brickLink link() const;
};

}

}

#endif
