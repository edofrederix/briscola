#ifndef decomposition_H
#define decomposition_H

#include "List.H"
#include "labelVector.H"
#include "block.H"
#include "brickDecompositionInterface.H"

namespace Foam
{

namespace briscola
{

class mesh;

class decomposition
{
protected:

    // Reference to the mesh
    mesh& msh_;

    // Reference to the decomposition dictionary
    const dictionary& dict_;

    // Decomposition per brick
    List<labelVector> decompPerBrick_;

    // The number of the brick of which this processor contains (part of) the
    // mesh
    label myBrickNum_;

    // The index which points to the part of this processor's brick which is
    // contained by this processor
    labelVector myBrickPart_;

private:

    // List of brick numbers per processor
    List<label> brickNumPerProc_;

    // List of part indices per processor
    List<labelVector> brickPartPerProc_;

    // List of part sizes per processor
    List<labelVector> partSizePerProc_;

    // List of processor topologies per brick
    PtrList<labelBlock> procMapPerBrick_;

    // List of part sizes per brick
    List<labelVector> partSizePerBrick_;

    // Global processor topology, only set on structured brick topologies
    labelBlock procMap_;

    // Copy constructor
    decomposition(const decomposition&);

    // Disallow default bitwise assignment
    void operator=(const decomposition&);

protected:

    // Update global data
    void updateGlobalData();

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

    // Selector
    static autoPtr<decomposition> New(mesh&);

    // Destructor
    virtual ~decomposition();

    // Decomposition dictionary
    inline const dictionary& dict() const
    {
        return dict_;
    }

    // Number of parts per brick
    inline const List<labelVector>& decompPerBrick() const
    {
        return decompPerBrick_;
    }

    // This processor's brick number
    inline label myBrickNum() const
    {
        return myBrickNum_;
    }

    // This processor's brick's shape
    labelVector myBrickN() const;

    // This processor's part's shape
    labelVector myPartN() const
    {
        return partSizePerBrick_[myBrickNum_];
    }

    // This processor's part index of its brick
    inline labelVector myBrickPart() const
    {
        return myBrickPart_;
    }

    // This processor's part start index within its brick
    inline labelVector myBrickPartStart() const
    {
        return cmptMultiply(myBrickPart_, partSizePerBrick_[myBrickNum_]);
    }

    // This processor's brick's number of parts
    inline labelVector myBrickDecomp() const
    {
        return decompPerBrick_[myBrickNum_];
    }

    // This processor's brick
    const brick& myBrick() const;

    // Brick numbers per processor
    inline const List<label>& brickNumPerProc() const
    {
        return brickNumPerProc_;
    }

    // Brick part per processor
    inline const List<labelVector>& brickPartPerProc() const
    {
        return brickPartPerProc_;
    }

    // Part sizes per processor
    inline const List<labelVector>& partSizePerProc() const
    {
        return partSizePerProc_;
    }

    // Global processor topology per brick
    inline const PtrList<labelBlock>& procMapPerBrick() const
    {
        return procMapPerBrick_;
    }

    // Global processor topology, only set on structured brick topologies
    inline const labelBlock& procMap() const
    {
        return procMap_;
    }

    // Global list of part sizes per brick
    inline const List<labelVector>& partSizePerBrick() const
    {
        return partSizePerBrick_;
    }

    // Given a part index and brick number, return the processor number
    inline label getProcNum(const labelVector p, const label b) const
    {
        forAll(brickPartPerProc_, proci)
        {
            if
            (
                brickPartPerProc_[proci] == p
             && brickNumPerProc_[proci] == b
            )
            {
                return proci;
            }
        }

        return -1;
    }
};

}

}

#endif
