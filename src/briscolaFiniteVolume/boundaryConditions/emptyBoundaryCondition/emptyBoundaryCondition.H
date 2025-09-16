#ifndef emptyBoundaryCondition_H
#define emptyBoundaryCondition_H

#include "emptyBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Forward declaration

template<class Type, class MeshType>
class emptyBoundaryCondition;


// Colocated

template<class Type>
class emptyBoundaryCondition<Type,colocated>
:
    public emptyBoundaryConditionBase<Type,colocated>
{
public:

    // Runtime type information
    TypeName("empty");

    // Use constructors from base
    using emptyBoundaryConditionBase<Type,colocated>::
        emptyBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,colocated>> clone() const
    {
        return tmp<boundaryCondition<Type,colocated>>
        (
            new emptyBoundaryCondition<Type,colocated>(*this)
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
            new emptyBoundaryCondition<Type,colocated>(field,*this)
        );
    }
};


// Staggered

template<class Type>
class emptyBoundaryCondition<Type,staggered>
:
    public emptyBoundaryConditionBase<Type,staggered>
{
public:

    // Runtime type information
    TypeName("empty");

    // Use constructors from base
    using emptyBoundaryConditionBase<Type,staggered>::
        emptyBoundaryConditionBase;

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new emptyBoundaryCondition<Type,staggered>(*this)
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
            new emptyBoundaryCondition<Type,staggered>(field,*this)
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
    #include "emptyBoundaryCondition.C"
#endif

#endif
