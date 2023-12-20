#ifndef stagNeumannBoundaryCondition_H
#define stagNeumannBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagNeumannBoundaryCondition
:
    public NeumannBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("Neumann");

    // Construct from components
    stagNeumannBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Construct from components and boundary gradient
    stagNeumannBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&,
        const List<Type>&
    );

    // Copy constructor
    stagNeumannBoundaryCondition
    (
        const stagNeumannBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagNeumannBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagNeumannBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagNeumannBoundaryCondition<Type>(*this)
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
            new stagNeumannBoundaryCondition<Type>(field,*this)
        );
    }

    // Eliminate ghosts from the linear system at given level and direction

    virtual void eliminateGhosts
    (
        linearSystem<stencil,Type,staggered>&,
        const label,
        const label
    );

    virtual void eliminateGhosts
    (
        linearSystem<symmStencil,Type,staggered>&,
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
    #include "stagNeumannBoundaryCondition.C"
#endif

#endif
