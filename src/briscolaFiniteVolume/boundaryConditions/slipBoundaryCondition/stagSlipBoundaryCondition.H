#ifndef stagSlipBoundaryCondition_H
#define stagSlipBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagSlipBoundaryCondition
:
    public NeumannBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("slip")

    // Construct from components
    stagSlipBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Copy constructor
    stagSlipBoundaryCondition
    (
        const stagSlipBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagSlipBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagSlipBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagSlipBoundaryCondition<Type>(*this)
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
            new stagSlipBoundaryCondition<Type>(field,*this)
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
    #include "stagSlipBoundaryCondition.C"
#endif

#endif
