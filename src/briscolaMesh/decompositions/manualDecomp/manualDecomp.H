#ifndef manualDecomp_H
#define manualDecomp_H

#include "decomposition.H"

namespace Foam
{

namespace briscola
{

namespace decompositions
{

/**

@brief Manual decomposition

This derived decomposition class requires a manual brick decomposition to be
provided by the user. For each brick the decomposition is to be specified in the
mesh dictionary. The processor distribution across bricks is then automatically
chosen by traversing the brick parts.

*/

class manualDecomp
:
    public decomposition
{
    // This processor's brick number
    label myBrickNum_;

    // This processor's brick decomposition
    labelVector myBrickDecomp_;

    // This processor's brick part
    labelVector myBrickPart_;

public:

    // Runtime type information
    TypeName("manual");

    // Construct from mesh
    manualDecomp(mesh&);

    // Copy constructor
    manualDecomp(const manualDecomp&);

    // Clone decomposition
    virtual autoPtr<decomposition> clone() const
    {
        return autoPtr<decomposition>(new manualDecomp(*this));
    }

    // This processor's brick number
    virtual label myBrickNum() const
    {
        return myBrickNum_;
    }

    // This processor's brick decomposition
    virtual labelVector myBrickDecomp() const
    {
        return myBrickDecomp_;
    }

    // This processor's part index of its brick
    virtual labelVector myBrickPart() const
    {
        return myBrickPart_;
    }
};

}

}

}

#endif
