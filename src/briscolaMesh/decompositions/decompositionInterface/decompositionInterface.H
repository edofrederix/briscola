#ifndef decompositionInterface_H
#define decompositionInterface_H

#include "brickLink.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "decompositionSlice.H"

namespace Foam
{

namespace briscola
{

class level;

class decompositionInterface
{
    // Copy constructor
    decompositionInterface(const decompositionInterface&);

    // Disallow default bitwise assignment
    void operator=(const decompositionInterface&);

protected:

    // Reference to the level
    const level& lvl_;

    // Brick link
    brickLink link_;

    // Reference to the first brick
    const brick& b0_;

    // Reference to the second brick
    const brick& b1_;

    // List of processor interface slices accross this brick interface
    PtrList<decompositionSlice> slices_;

    // Slice topology
    labelBlock map_;

    // Set the offset vector
    void setOffsetAndTransformation();

    // Set the processor interface slices
    void setSlices();

public:

    // Construct from brick link and level
    decompositionInterface(const brickLink&, const level&);

    // Destructor
    virtual ~decompositionInterface();

    // Return access the brick link
    inline const brickLink& link() const { return link_; }

    // Return access to the list of slices
    const PtrList<decompositionSlice>& slices() const { return slices_; }

    // Return access to the slice topology
    const labelBlock& map() const { return map_; }
};

}

}

#endif
