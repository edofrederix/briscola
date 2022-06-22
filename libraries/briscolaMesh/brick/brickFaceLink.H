#ifndef brickFaceLink_H
#define brickFaceLink_H

#include "face.H"
#include "labelVector.H"
#include "labelTensor.H"

namespace Foam
{

namespace briscola
{

class brickTopology;

class brickFaceLink
{
    // Reference to the parent face
    const face& f0_;

    // Reference to the daughter face
    const face& f1_;

    // Transformation tensor of the daughter brick to the parent brick
    labelTensor T_;

    // Periodic link
    bool periodic_;

public:

    // Construct from parent and daughter Face. Optionally set as periodic link
    brickFaceLink(const face&, const face&, const bool = false);

    // Construct from parent, daughter, path and topology. Optionally set as
    // periodic link
    brickFaceLink
    (
        const face&,
        const face&,
        const labelList&,
        const brickTopology&,
        const bool = false
    );

    // Destructor
    virtual ~brickFaceLink();

    // Access to the parent face
    inline const face& f0() const { return f0_; }

    // Access to the daughter face
    inline const face& f1() const { return f1_; }

    // Calculate the offset
    inline labelVector offset() const
    {
        return faceOffsets[f0_.num()];
    }

    // Access to the transformation tensor
    inline const labelTensor& T() const { return T_; }

    // Periodic link
    inline bool periodic() const { return periodic_; }
};

}

}

#endif
