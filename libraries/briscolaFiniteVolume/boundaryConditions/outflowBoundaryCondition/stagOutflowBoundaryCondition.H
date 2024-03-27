#ifndef stagOutflowBoundaryCondition_H
#define stagOutflowBoundaryCondition_H

#include "outflowBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagOutflowBoundaryCondition
:
    public outflowBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("outflow");

    // Construct from components
    stagOutflowBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Copy constructor
    stagOutflowBoundaryCondition
    (
        const stagOutflowBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagOutflowBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagOutflowBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagOutflowBoundaryCondition<Type>(*this)
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
            new stagOutflowBoundaryCondition<Type>(field,*this)
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
    #include "stagOutflowBoundaryCondition.C"
#endif

#endif
