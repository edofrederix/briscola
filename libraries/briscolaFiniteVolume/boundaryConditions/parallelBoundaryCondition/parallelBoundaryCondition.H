#ifndef parallelBoundaryCondition_H
#define parallelBoundaryCondition_H

#include "boundaryCondition.H"

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
    PtrList<block<Type>> sendBuffers_;

    // Receive buffers for each mesh direction
    PtrList<block<Type>> recvBuffers_;

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
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    parallelBoundaryCondition(const parallelBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    parallelBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const parallelBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new parallelBoundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new parallelBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return PARALLELBC;
    }

    // Non-eliminated boundary condition
    virtual bool eliminated() const { return false; }

    // Prepare the evaluation of the boundary condition for given level
    virtual void prepare(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);
};

}

}

}

#ifdef NoRepository
    #include "parallelBoundaryCondition.C"
#endif

#endif
