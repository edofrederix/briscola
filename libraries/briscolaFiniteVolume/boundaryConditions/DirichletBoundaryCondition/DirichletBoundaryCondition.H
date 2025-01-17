#ifndef DirichletBoundaryCondition_H
#define DirichletBoundaryCondition_H

#include "boundaryCondition.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class DirichletBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
protected:

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("Dirichlet");

    // Construct from components
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and boundary values
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
    );

    // Copy constructor
    DirichletBoundaryCondition
    (
        const DirichletBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const DirichletBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new DirichletBoundaryCondition<Type,MeshType>(*this)
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
            new DirichletBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DIRICHLETBC;
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

    // Evaluate the boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "DirichletBoundaryCondition.C"
#endif

#endif
