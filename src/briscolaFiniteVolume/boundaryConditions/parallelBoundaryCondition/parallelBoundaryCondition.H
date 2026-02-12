#ifndef parallelBoundaryCondition_H
#define parallelBoundaryCondition_H

#include "boundaryCondition.H"
#include "parallelBoundary.H"

#include "block.H"
#include "faceLabel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class parallelBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
protected:

    // Neighbor processor number
    label neighborProcNum_;

    // Communication tag
    label tag_;

    // Send buffers for each mesh direction
    FastPtrList<block<Type>> sendBuffers_;

    // Receive buffers for each mesh direction
    FastPtrList<block<Type>> recvBuffers_;

    // Outstanding send requests
    label outstandingSendRequest_;

    // Outstanding receive requests
    label outstandingRecvRequest_;

private:

    // Disallow default bitwise assignment
    void operator=(const parallelBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("parallel");

    // Construct from components
    parallelBoundaryCondition
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    parallelBoundaryCondition(const parallelBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting level reference
    parallelBoundaryCondition
    (
        const parallelBoundaryCondition<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new parallelBoundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshLevel<Type,MeshType>& level
    ) const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new parallelBoundaryCondition<Type,MeshType>(*this, level)
        );
    }

    // Neighbor processor number
    label neighborProcNum() const { return neighborProcNum_; }

    // Tag number
    label tag() const { return tag_; }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return PARALLELBC;
    }

    // Return the transformation tensor
    inline labelTensor T() const
    {
        return this->b_.template cast<parallelBoundary>().T();
    }

    // Non-eliminated boundary condition
    virtual bool eliminated() const { return false; }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();

    // Evaluate the boundary condition
    virtual void evaluate();
};

}

}

}

#ifdef NoRepository
    #include "parallelBoundaryCondition.C"
#endif

#endif
