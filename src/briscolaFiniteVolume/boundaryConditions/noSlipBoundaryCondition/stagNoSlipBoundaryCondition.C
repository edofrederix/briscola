#include "stagNoSlipBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagNoSlipBoundaryCondition<Type>::stagNoSlipBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    DirichletBoundaryCondition<Type,staggered>
    (
        mshField,
        b,
        List<Type>(staggered::numberOfDirections, Zero)
    )
{}

template<class Type>
stagNoSlipBoundaryCondition<Type>::stagNoSlipBoundaryCondition
(
    const stagNoSlipBoundaryCondition<Type>& bc
)
:
    DirichletBoundaryCondition<Type,staggered>(bc)
{}

template<class Type>
stagNoSlipBoundaryCondition<Type>::stagNoSlipBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagNoSlipBoundaryCondition<Type>& bc
)
:
    DirichletBoundaryCondition<Type,staggered>(field, bc)
{}

}

}

}

