#ifndef periodicBoundaryCondition_H
#define periodicBoundaryCondition_H

#include "parallelBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class periodicBoundaryCondition
:
    public parallelBoundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicBoundaryCondition
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    periodicBoundaryCondition(const periodicBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting level reference
    periodicBoundaryCondition
    (
        const periodicBoundaryCondition<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new periodicBoundaryCondition<Type,MeshType>(*this)
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
            new periodicBoundaryCondition<Type,MeshType>(*this, level)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return PERIODICBC;
    }

    // Prepare the evaluation of the boundary condition
    virtual void prepare();

    // Evaluate the boundary condition
    virtual void evaluate();
};

}

}

}

#ifdef NoRepository
    #include "periodicBoundaryCondition.C"
#endif

#endif
