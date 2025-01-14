#ifndef brickDecompositionInterface_H
#define brickDecompositionInterface_H

#include "brickLink.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "brickDecompositionSlice.H"

namespace Foam
{

namespace briscola
{

class mesh;

class brickDecompositionInterface
{
    // Copy constructor
    brickDecompositionInterface(const brickDecompositionInterface&);

    // Disallow default bitwise assignment
    void operator=(const brickDecompositionInterface&);

protected:

    // Brick link
    brickLink link_;

    // Reference to the first brick
    const brick& b0_;

    // Reference to the second brick
    const brick& b1_;

    // Reference to the mesh
    const mesh& msh_;

    // List of processor interface slices accross this brick interface
    PtrList<brickDecompositionSlice> slices_;

    // Slice topology
    labelBlock map_;

    // Set the offset vector
    void setOffsetAndTransformation();

    // Set the processor interface slices
    void setSlices();

public:

    // Construct from brick link and mesh
    brickDecompositionInterface(const brickLink&, const mesh&);

    // Destructor
    virtual ~brickDecompositionInterface();

    // Return access the brick link
    inline const brickLink& link() const { return link_; }

    // Return access to the list of slices
    const PtrList<brickDecompositionSlice>& slices() const { return slices_; }

    // Return access to the slice topology
    const labelBlock& map() const { return map_; }
};

}

}

#endif
