#ifndef stagEmptyBoundaryCondition_H
#define stagEmptyBoundaryCondition_H

#include "emptyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagEmptyBoundaryCondition
:
    public emptyBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("empty");

    // Construct from components
    stagEmptyBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Copy constructor
    stagEmptyBoundaryCondition
    (
        const stagEmptyBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagEmptyBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagEmptyBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagEmptyBoundaryCondition<Type>(*this)
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
            new stagEmptyBoundaryCondition<Type>(field,*this)
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
    #include "stagEmptyBoundaryCondition.C"
#endif

#endif
