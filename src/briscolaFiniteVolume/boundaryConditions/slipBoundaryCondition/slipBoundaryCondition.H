#ifndef slipBoundaryCondition_H
#define slipBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class slipBoundaryCondition;


// Colocated

template<class Type>
class slipBoundaryCondition<Type,colocated>
:
    public boundaryCondition<Type,colocated>
{
public:

    // Runtime type information
    TypeName("slip");

    // Use constructors from base
    using boundaryCondition<Type,colocated>::boundaryCondition;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new slipBoundaryCondition<Type,colocated>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<Type,colocated>> clone
    (
        const meshLevel<Type,colocated>& level
    ) const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new slipBoundaryCondition<Type,colocated>(*this, level)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return false; }

    // Evaluate the boundary condition for given direction
    virtual void evaluate(const label);
};


// Staggered

template<class Type>
class slipBoundaryCondition<Type,staggered>
:
    public NeumannBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("slip")

    // Use constructors from base
    using NeumannBoundaryCondition<Type,staggered>::
        NeumannBoundaryCondition;

    // Construct from components setting zero gradient
    slipBoundaryCondition
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
            new slipBoundaryCondition<Type,staggered>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<Type,staggered>> clone
    (
        const meshLevel<Type,staggered>& level
    ) const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new slipBoundaryCondition<Type,staggered>(*this, level)
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
    #include "slipBoundaryCondition.C"
#endif

#endif
