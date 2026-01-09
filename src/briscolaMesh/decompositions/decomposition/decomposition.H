#ifndef decomposition_H
#define decomposition_H

#include "List.H"
#include "FastPtrList.H"
#include "labelVector.H"
#include "block.H"
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

/**

@brief Decomposition class

This is the base class defining the decomposition of the level. This base class
handles general collection of data. The actual decomposition must be defined by
a derived type, by means of the myBrickNum(), myBrickDecomp() and myBrickPart()
pure virtual functions. On each processor, these functions return the brick
number to which the processor belongs, the decomposition of that brick and the
brick part of that brick to which the processor belongs. Based on this
information, the decomposition class fills some general data arrays. This is
done in the init() function which is automatically called in the selector. If
the selector is not used, then the init() function must be called manually.

*/

class level;

class decomposition
{
protected:

    // Reference to the level
    const level& lvl_;

    // Reference to the decomposition dictionary
    const dictionary& dict_;

private:

    // Data for each processor

        // Processor brick numbers
        labelList procBrickNums_;

        // Processor brick decompositions
        List<labelVector> procBrickDecomps_;

        // Processor brick parts
        List<labelVector> procBrickParts_;

    // Data for each brick

        // Brick decompositions
        List<labelVector> brickDecomps_;

        // Brick processor topologies
        FastPtrList<labelBlock> brickProcMaps_;

    // Structured processor topology map (only set on structured meshes)
    autoPtr<decompositionMap> mapPtr_;

public:

    // Runtime type information
    TypeName("decomposition");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        decomposition,
        dictionary,
        (const level& lvl),
        (lvl)
    );

    // Construct from level
    decomposition(const level&);

    // Copy constructor
    decomposition(const decomposition&);

    // Clone decomposition
    virtual autoPtr<decomposition> clone() const = 0;

    // Selector
    static autoPtr<decomposition> New(const level&);

    // Destructor
    virtual ~decomposition();

    // Level
    inline const level& lvl() const
    {
        return lvl_;
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

    // Initialize
    void init();

    // Pure virtual functions to be defined by derived types

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
