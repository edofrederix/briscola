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
        const meshLevel<Type,colocated>& level,
        const boundary& b
    )
    :
        NeumannBoundaryCondition<Type,colocated>
        (
            level,
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
        const meshLevel<Type,staggered>& level,
        const boundary& b
    )
    :
        NeumannBoundaryCondition<Type,staggered>
        (
            level,
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
    #include "outflowBoundaryCondition.C"
#endif

#endif
