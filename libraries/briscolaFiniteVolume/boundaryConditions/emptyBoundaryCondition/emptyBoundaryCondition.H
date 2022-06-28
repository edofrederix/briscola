#ifndef emptyBoundaryCondition_H
#define emptyBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
private:

    // Disallow default bitwise assignment
    void operator=(const emptyBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("empty");

    // Construct from components
    emptyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Copy constructor
    emptyBoundaryCondition(const emptyBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    emptyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const emptyBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new emptyBoundaryCondition<Type,MeshType>(*this)
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
            new emptyBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return EMPTYBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label, const bool = false);

    // Empty boundary condition are non-eliminated on shifted boundaries
    virtual bool eliminated() const { return false; }
};

}

}

}

#endif
