#ifndef stagDirichletBoundaryCondition_H
#define stagDirichletBoundaryCondition_H

#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagDirichletBoundaryCondition
:
    public DirichletBoundaryCondition<Type,staggered>
{
public:

    // Runtime type information
    TypeName("Dirichlet");

    // Construct from components
    stagDirichletBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Construct from components and boundary values
    stagDirichletBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&,
        const List<Type>&
    );

    // Copy constructor
    stagDirichletBoundaryCondition
    (
        const stagDirichletBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagDirichletBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagDirichletBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagDirichletBoundaryCondition<Type>(*this)
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
            new stagDirichletBoundaryCondition<Type>(field,*this)
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
    #include "stagDirichletBoundaryCondition.C"
#endif

#endif
