#ifndef outflowBoundaryCondition_H
#define outflowBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class outflowBoundaryCondition;


// Colocated

template<class Type>
class outflowBoundaryCondition<Type,colocated>
:
    public NeumannBoundaryCondition<Type,colocated>
{
public:

    // Runtime type information
    TypeName("outflow");

    // Use base class constructors
    using NeumannBoundaryCondition<Type,colocated>::
        NeumannBoundaryCondition;

    // Construct from components
    outflowBoundaryCondition
    (
        const meshField<Type,colocated>& mshField,
        const boundary& b
    )
    :
        NeumannBoundaryCondition<Type,colocated>
        (
            mshField,
            b,
            Zero
        )
    {}

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new outflowBoundaryCondition<Type,colocated>(*this)
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
            new outflowBoundaryCondition<Type,colocated>(field,*this)
        );
    }
};


// Staggered

template<class Type>
class outflowBoundaryCondition<Type,staggered>
:
    public NeumannBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("outflow");

    // Use base class constructors
    using NeumannBoundaryCondition<Type,staggered>::
        NeumannBoundaryCondition;

    // Construct from components
    outflowBoundaryCondition
    (
        const meshField<Type,staggered>& mshField,
        const boundary& b
    )
    :
        NeumannBoundaryCondition<Type,staggered>
        (
            mshField,
            b,
            Zero
        )
    {}

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new outflowBoundaryCondition<Type,staggered>(*this)
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
            new outflowBoundaryCondition<Type,staggered>(field,*this)
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
    #include "outflowBoundaryCondition.C"
#endif

#endif
