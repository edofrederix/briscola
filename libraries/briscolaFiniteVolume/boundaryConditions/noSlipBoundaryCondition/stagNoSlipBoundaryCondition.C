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
    stagDirichletBoundaryCondition<Type>
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
    stagDirichletBoundaryCondition<Type>
    (
        bc.mshField(),
        bc.mshBoundary()
    )
{}

template<class Type>
stagNoSlipBoundaryCondition<Type>::stagNoSlipBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagNoSlipBoundaryCondition<Type>& bc
)
:
    stagDirichletBoundaryCondition<Type>
    (
        field,
        bc.mshBoundary()
    )
{}

}

}

}

