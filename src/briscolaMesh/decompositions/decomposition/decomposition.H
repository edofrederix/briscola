#ifndef decomposition_H
#define decomposition_H

#include "List.H"
#include "FastPtrList.H"
#include "labelVector.H"
#include "block.H"
#include "brickDecompositionInterface.H"
#include "decompositionMap.H"

#include "faceLabel.H"
#include "faceLabelTensor.H"
#include "edgeLabel.H"
#include "edgeLabelTensor.H"
#include "vertexLabel.H"
#include "vertexLabelTensor.H"

namespace Foam
{

namespace briscola
{

class mesh;

class decomposition
{
protected:

    // Reference to the mesh
    const mesh& msh_;

    // Reference to the decomposition dictionary
    const dictionary& dict_;

private:

    // Processor brick numbers
    labelList procBrickNums_;

    // Processor brick decompositions
    List<labelVector> procBrickDecomps_;

    // Processor brick parts
    List<labelVector> procBrickParts_;

    // Brick decompositions
    List<labelVector> brickDecomps_;

    // Brick processor topologies
    FastPtrList<labelBlock> brickProcMaps_;

    // Global processor topology map, only set on structured brick topologies
    autoPtr<decompositionMap> mapPtr_;

protected:

    // Initialize
    void init();

public:

    // Runtime type information
    TypeName("decomposition");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        decomposition,
        dictionary,
        (
            mesh& msh
        ),
        (msh)
    );

    // Construct from mesh
    decomposition(mesh&);

    // Copy constructor
    decomposition(const decomposition&);

    // Clone decomposition
    virtual autoPtr<decomposition> clone() const = 0;

    // Selector
    static autoPtr<decomposition> New(mesh&);

    // Destructor
    virtual ~decomposition();

    // Mesh
    inline const mesh& msh() const
    {
        return msh_;
    }

    // Decomposition dictionary
    inline const dictionary& dict() const
    {
        return dict_;
    }

    // Processor brick numbers
    inline const labelList& procBrickNums() const
    {
        return procBrickNums_;
    }

    // Processor brick decompositions
    inline const List<labelVector>& procBrickDecomps() const
    {
        return procBrickDecomps_;
    }

    // Processor brick parts
    inline const List<labelVector>& procBrickParts() const
    {
        return procBrickParts_;
    }

    // Brick decompositions
    inline const List<labelVector>& brickDecomps() const
    {
        return brickDecomps_;
    }

    // Brick processor topologies
    inline const FastPtrList<labelBlock>& brickProcMaps() const
    {
        return brickProcMaps_;
    }

    // Global processor topology, only set on structured brick topologies
    inline const decompositionMap& map() const
    {
        return mapPtr_();
    }

    // This processor's brick number
    virtual label myBrickNum() const = 0;

    // This processor's brick decomposition
    virtual labelVector myBrickDecomp() const = 0;

    // This processor's part index of its brick
    virtual labelVector myBrickPart() const = 0;
};

}

}

#endif
