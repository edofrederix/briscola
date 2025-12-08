#ifndef decompositionMap_H
#define decompositionMap_H

#include "labelBlock.H"
#include "List.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

class decomposition;

class decompositionMap
:
    public labelBlock
{
    // Reference to the decomposition
    const decomposition& decomp_;

    // Brick legend
    List<labelVector> legend_;

    // Part shapes
    labelVectorBlock shapes_;

    // Part starts
    labelVectorBlock starts_;

public:

    // Construct from decomposition
    decompositionMap(const decomposition&);

    // Copy constructor
    decompositionMap(const decompositionMap&);

    // Disallow null constructor
    decompositionMap() = delete;

    // Destructor
    virtual ~decompositionMap();

    // Access to the brick legend
    const List<labelVector>& legend() const { return legend_; }

    // Part shapes
    const labelVectorBlock& shapes() const { return shapes_; }

    // Part starts
    const labelVectorBlock& starts() const { return starts_; }

    // This processor's part shape
    labelVector myBrickPartShape() const
    {
        return shapes_(legend_[Pstream::myProcNo()]);
    }

    // This processor's part start
    labelVector myBrickPartStart() const
    {
        return starts_(legend_[Pstream::myProcNo()]);
    }
};

}

}

#endif
