#ifndef parallelBoundaryCondition_H
#define parallelBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

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

    // Order in which the mesh directions are received. Given the neighbor's
    // direction, order will return the local corresponding direction
    labelList order_;

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
        const partPatch&
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

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label, const bool = false);

    // Parallel face boundary conditions are non-eliminated on shifted slave
    // boundaries
    virtual bool eliminated() const
    {
        return !(this->slave() && this->boundaryOffsetDegree() == 1);
    }
};

}

}

}

#endif
