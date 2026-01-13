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

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,colocated>> clone
    (
        const meshField<Type,colocated>& field
    ) const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new slipBoundaryCondition<Type,colocated>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return false; }

    // Evaluate the boundary condition for given level and direction
    virtual void evaluate(const label, const label);
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
            new slipBoundaryCondition<Type,staggered>(*this)
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
            new slipBoundaryCondition<Type,staggered>(field,*this)
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
    #include "slipBoundaryCondition.C"
#endif

#endif
