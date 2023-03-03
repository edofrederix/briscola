#ifndef brickEdgeLink_H
#define brickEdgeLink_H

#include "edge.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "labelList.H"

namespace Foam
{

namespace briscola
{

class brickTopology;

class brickEdgeLink
{
    // Reference to the parent edge
    const edge& e0_;

    // Reference to the daughter edge
    const edge& e1_;

    // Transformation tensor of the daughter brick to the parent brick
    labelTensor T_;

    // Periodic link
    bool periodic_;

public:

    // Construct from parent, daughter and transformation tensor. Optionally set
    // as periodic link.
    brickEdgeLink
    (
        const edge&,
        const edge&,
        const labelTensor&,
        const bool = false
    );

    // Copy constructor
    inline brickEdgeLink(const brickEdgeLink& link)
    :
        e0_(link.e0_),
        e1_(link.e1_),
        T_(link.T_),
        periodic_(link.periodic_)
    {}

    // Clone link
    inline autoPtr<brickEdgeLink> clone() const
    {
        return autoPtr<brickEdgeLink>(new brickEdgeLink(*this));
    }

    // Destructor
    virtual ~brickEdgeLink();

    // Access to the parent edge
    inline const edge& e0() const { return e0_; }

    // Access to the daughter edge
    inline const edge& e1() const { return e1_; }

    // Calculate the offset
    inline labelVector offset() const
    {
        return edgeOffsets[e0_.num()];
    }

    // Access to the transformation tensor
    inline const labelTensor& T() const { return T_; }

    // Periodic link
    inline bool periodic() const { return periodic_; }
};

}

}

#endif
