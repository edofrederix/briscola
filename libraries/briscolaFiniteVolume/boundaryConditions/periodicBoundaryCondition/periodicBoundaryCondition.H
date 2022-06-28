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
    // Disallow default bitwise assignment
    void operator=(const periodicBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("periodic");

    // Construct from components
    periodicBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Copy constructor
    periodicBoundaryCondition(const periodicBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    periodicBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const periodicBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new periodicBoundaryCondition<Type,MeshType>(*this)
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
            new periodicBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return PERIODICBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);
};

}

}

}

#endif
