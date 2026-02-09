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

    // Agglomerate data (if this is an agglomerate member)

        // Our agglomerate processor map
        labelBlock aggProcMap_;

        // Our agglomerate part
        labelVector myAggPart_;

        // Our agglomerate processor number
        label myAggProcNo_;

    // Members of this decomposition
    labelList members_;

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

    // Copy constructor with alternative level reference
    decomposition(const decomposition&, const level&);

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

    // Our agglomerate map
    const labelBlock& aggProcMap() const
    {
        return aggProcMap_;
    }

    // Our agglomerate part
    labelVector myAggPart() const
    {
        return myAggPart_;
    }

    // Our agglomerate processor number (local value)
    label myAggProcNo() const
    {
        return myAggProcNo_;
    }

    // Global processor topology, only set on structured brick topologies
    inline const decompositionMap& map() const
    {
        return mapPtr_();
    }

    // Initialize
    void init();

    // Check if decomposition is coarsenable, which is true if we have a
    // non-unity decomposition and if the agglomerate mesh can be coarsened
    bool coarsenable() const;

    // Pure virtual functions to be defined by derived types

        // This processor's brick number
        virtual label myBrickNum() const = 0;

        // This processor's brick decomposition
        virtual labelVector myBrickDecomp() const = 0;

        // This processor's part index of its brick
        virtual labelVector myBrickPart() const = 0;

    // Function to check whether this processor is member of the decomposition
    inline bool member() const
    {
        return myBrickPart() != -unitXYZ;
    }

    // Return list of member processors
    inline const labelList& members() const { return members_; }

    // Agglomeration ratio
    inline labelVector R() const
    {
        return cmptMax(aggProcMap_.shape(), unitXYZ);
    }

    // A decomposition is agglomerated if there is a non-unity agglomeration
    // ratio and if we are a member of the parent decomposition
    bool agglomerated() const;

    // A level is an agglomerate master if it is agglomerated and a member of
    // this level's decomposition
    inline bool aggMaster() const
    {
        return agglomerated() && member();
    }

    // A level is an agglomerate slave if it is agglomerated and not a member of
    // this level's decomposition
    inline bool aggSlave() const
    {
        return agglomerated() && !member();
    }

    // A level is an agglomerate parent if its child decomposition is either an
    // agglomerate master of slave
    bool aggParent() const;
};

}

}

#endif
