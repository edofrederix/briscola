#ifndef DirichletBoundaryCondition_H
#define DirichletBoundaryCondition_H

#include "DirichletBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class DirichletBoundaryCondition;


// Colocated

template<class Type>
class DirichletBoundaryCondition<Type,colocated>
:
    public DirichletBoundaryConditionBase<Type,colocated>
{
public:

    // Runtime type information
    TypeName("Dirichlet");

    // Use constructors from base
    using DirichletBoundaryConditionBase<Type,colocated>::
        DirichletBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new DirichletBoundaryCondition<Type,colocated>(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,colocated>> clone
    (
        const meshField<Type,colocated>& field
    ) const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new DirichletBoundaryCondition<Type,colocated>(field,*this)
        );
    }

    // Eliminate ghosts from the linear system at given level and direction
    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,colocated>&,
        const label,
        const label
    );

    // Evaluate the boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};


// Staggered

template<class Type>
class DirichletBoundaryCondition<Type,staggered>
:
    public DirichletBoundaryConditionBase<Type,staggered>
{
public:

    // Runtime type information
    TypeName("Dirichlet");

    // Use constructors from base
    using DirichletBoundaryConditionBase<Type,staggered>::
        DirichletBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new DirichletBoundaryCondition<Type,staggered>(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,staggered>> clone
    (
        const meshField<Type,staggered>& field
    ) const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new DirichletBoundaryCondition<Type,staggered>(field,*this)
        );
    }

    // Eliminate ghosts from the linear system at given level and direction
    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,staggered>&,
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
