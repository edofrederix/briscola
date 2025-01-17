#ifndef brickLink_H
#define brickLink_H

#include "labelVector.H"
#include "labelTensor.H"
#include "brick.H"

#include "brickFaceLink.H"
#include "brickEdgeLink.H"
#include "brickVertexLink.H"

namespace Foam
{

namespace briscola
{

class brickLink
{
    // Reference to the parent brick
    const brick& b0_;

    // Reference to the daughter brick
    const brick& b1_;

    // Offset vector from the parent to the daughter brick
    const labelVector offset_;

    // Transformation tensor of the daughter brick to the parent brick
    const labelTensor T_;

    // Periodic link
    bool periodic_;

public:

    // Construct from components
    brickLink
    (
        const brick&,
        const brick&,
        const labelVector = zeroXYZ,
        const labelTensor = eye,
        const bool = false
    );

    // Copy constructor
    brickLink(const brickLink&);

    // Construct from face link
    brickLink(const brickFaceLink&);

    // Construct from edge link
    brickLink(const brickEdgeLink&);

    // Construct from vertex link
    brickLink(const brickVertexLink&);

    // Destructor
    virtual ~brickLink();

    // Access to the parent brick
    inline const brick& b0() const { return b0_; }

    // Access to the daughter brick
    inline const brick& b1() const { return b1_; }

    // Return the offset vector
    inline labelVector offset() const { return offset_; }

    // Return the transformation tensor
    inline labelTensor T() const { return T_; }

    // Periodic
    inline bool periodic() const { return periodic_; }
};

}

}

#endif
