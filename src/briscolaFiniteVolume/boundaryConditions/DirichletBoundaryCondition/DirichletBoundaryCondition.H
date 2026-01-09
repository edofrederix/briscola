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

    // Eliminate ghosts from the linear system at given direction
    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,colocated>&,
        const label
    );

    // Evaluate the boundary condition for given direction
    virtual void evaluate(const label);
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

    // Eliminate ghosts from the linear system at given direction
    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,staggered>&,
        const label
    );

    // Evaluate the boundary condition for given direction
    virtual void evaluate(const label);
};

}

}

}

#ifdef NoRepository
    #include "DirichletBoundaryCondition.C"
#endif

#endif
