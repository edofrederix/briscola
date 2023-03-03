#ifndef decomposition_H
#define decomposition_H

#include "List.H"
#include "labelVector.H"
#include "block.H"
#include "brickDecompositionInterface.H"
#include "decompositionMap.H"

namespace Foam
{

namespace briscola
{

class mesh;

class decomposition
{
protected:

    // Reference to the decomposition dictionary
    const dictionary& dict_;

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

    // Global processor topology map, only set on structured brick topologies
    decompositionMap map_;

    // Disallow default bitwise assignment
    void operator=(const decomposition&);

protected:

    // Update global data
    void updateGlobalData(const mesh&);

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

    // Decomposition dictionary
    inline const dictionary& dict() const
    {
        return dict_;
    }

    // Number of parts per brick, to be determined by derived class
    virtual const List<labelVector>& decompPerBrick() const = 0;

    // This processor's brick number, to be determined by derived class
    virtual label myBrickNum() const = 0;

    // This processor's part index of its brick, to be determined by derived
    // class
    virtual labelVector myBrickPart() const = 0;

    // This processor's part's shape
    labelVector myPartN() const
    {
        return partSizePerBrick_[myBrickNum()];
    }

    // This processor's part start index within its brick
    inline labelVector myBrickPartStart() const
    {
        return cmptMultiply(myBrickPart(), partSizePerBrick_[myBrickNum()]);
    }

    // This processor's brick's number of parts
    inline labelVector myBrickDecomp() const
    {
        return decompPerBrick()[myBrickNum()];
    }

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
    inline const decompositionMap& map() const
    {
        return map_;
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
