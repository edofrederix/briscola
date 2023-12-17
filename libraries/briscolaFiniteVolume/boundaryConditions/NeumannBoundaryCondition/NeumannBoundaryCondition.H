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

    // Boundary gradients
    List<Type> boundaryGradients_;

public:

    // Runtime type information
    TypeName("Neumann");

    // Construct from components
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and boundary gradient
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
    );

    // Copy constructor
    NeumannBoundaryCondition
    (
        const NeumannBoundaryCondition<Type,MeshType>&
    );

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

    // Eliminated boundary condition
    virtual bool eliminated() const { return true; }

    // Eliminate ghosts from the linear system at given level and direction

    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,MeshType>&,
        const label,
        const label
    );

    virtual void eliminateGhosts
    (
        linearSystem<symmStencil,Type,MeshType>&,
        const label,
        const label
    );

    // Evaluate the boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "NeumannBoundaryCondition.C"
#endif

#endif
