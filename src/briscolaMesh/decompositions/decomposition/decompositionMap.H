#ifndef decompositionMap_H
#define decompositionMap_H

#include "labelBlock.H"
#include "List.H"
#include "labelVector.H"

namespace Foam
{

namespace briscola
{

/**

@brief Structured decomposition map

This class contains a label block that represents the structured processor
topology. Each entry in the block gives the processor number that corresponds to
the brick part. It is -1 if there is no processor associated to with the part,
which happens for meshes which do not have a rectangular brick structure. The
legend gives for each processor the decomposition map index that it belongs to.

*/

class decomposition;

class decompositionMap
:
    public labelBlock
{
    // Reference to the decomposition
    const decomposition& decomp_;

    // Decomposition map legend
    DynamicList<labelVector> legend_;

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

    // Access to the decomposition map legend
    const DynamicList<labelVector>& legend() const { return legend_; }

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
