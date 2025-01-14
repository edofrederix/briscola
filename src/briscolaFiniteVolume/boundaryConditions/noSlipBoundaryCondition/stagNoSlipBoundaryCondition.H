#ifndef stagNoSlipBoundaryCondition_H
#define stagNoSlipBoundaryCondition_H

#include "stagDirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class stagNoSlipBoundaryCondition
:
    public stagDirichletBoundaryCondition<Type>
{
public:

    // Runtime type information
    TypeName("noSlip");

    // Construct from components
    stagNoSlipBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const boundary&
    );

    // Copy constructor
    stagNoSlipBoundaryCondition
    (
        const stagNoSlipBoundaryCondition<Type>&
    );

    // Copy constructor resetting field reference
    stagNoSlipBoundaryCondition
    (
        const meshField<Type,staggered>&,
        const stagNoSlipBoundaryCondition<Type>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,staggered>> clone() const
    {
        return tmp<boundaryCondition<Type,staggered>>
        (
            new stagNoSlipBoundaryCondition<Type>(*this)
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
            new stagNoSlipBoundaryCondition<Type>(field,*this)
        );
    }
};

}

}

}

#ifdef NoRepository
    #include "stagNoSlipBoundaryCondition.C"
#endif

#endif
