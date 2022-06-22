#ifndef NeumannBoundaryCondition_H
#define NeumannBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class NeumannBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
protected:

    // Boundary gradients for each level and direction
    PtrList<block<Type>> boundaryGradients_;

private:

    // Disallow default bitwise assignment
    void operator=(const NeumannBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("Neumann");

    // Construct from components
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Copy constructor
    NeumannBoundaryCondition(const NeumannBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const NeumannBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new NeumannBoundaryCondition<Type,MeshType>(*this)
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
            new NeumannBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label, const bool = false);

    // Neumann boundary condition does not constrain on-boundary values
    virtual bool constrained() const { return false; }
};

}

}

}

#endif
