#ifndef manualDecomp_H
#define manualDecomp_H

#include "decomposition.H"

namespace Foam
{

namespace briscola
{

namespace decompositions
{

class manualDecomp
:
    public decomposition
{
    // Decomposition per brick
    List<labelVector> decompPerBrick_;

    // The number of the brick of which this processor contains (part of) the
    // mesh
    label myBrickNum_;

    // The index which points to the part of this processor's brick which is
    // contained by this processor
    labelVector myBrickPart_;

    // Copy constructor
    manualDecomp(const manualDecomp&);

    // Disallow default bitwise assignment
    void operator=(const manualDecomp&);

public:

    // Number of parts per brick
    virtual const List<labelVector>& decompPerBrick() const
    {
        return decompPerBrick_;
    }

    // This processor's brick number
    virtual label myBrickNum() const
    {
        return myBrickNum_;
    }

    // This processor's part index of its brick
    virtual labelVector myBrickPart() const
    {
        return myBrickPart_;
    }

    // Runtime type information
    TypeName("manual");

    // Construct from mesh
    manualDecomp(mesh&);

};

}

}

}

#endif
