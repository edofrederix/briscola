#ifndef NeumannBoundaryCondition_H
#define NeumannBoundaryCondition_H

#include "NeumannBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class NeumannBoundaryCondition;


// Colocated

template<class Type>
class NeumannBoundaryCondition<Type,colocated>
:
    public NeumannBoundaryConditionBase<Type,colocated>
{
public:

    // Runtime type information
    TypeName("Neumann");

    // Use constructors from base
    using NeumannBoundaryConditionBase<Type,colocated>::
        NeumannBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new NeumannBoundaryCondition<Type,colocated>(*this)
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
            new NeumannBoundaryCondition<Type,colocated>(field,*this)
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
class NeumannBoundaryCondition<Type,staggered>
:
    public NeumannBoundaryConditionBase<Type,staggered>
{
public:

    // Runtime type information
    TypeName("Neumann");

    // Use constructors from base
    using NeumannBoundaryConditionBase<Type,staggered>::
        NeumannBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new NeumannBoundaryCondition<Type,staggered>(*this)
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
            new NeumannBoundaryCondition<Type,staggered>(field,*this)
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
    #include "NeumannBoundaryCondition.C"
#endif

#endif
