#include "stagOutflowBoundaryCondition.H"

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
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const meshField<Type,staggered>& mshField,
    const boundary& b
)
:
    outflowBoundaryCondition<Type,staggered>(mshField, b)
{}

template<class Type>
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const stagOutflowBoundaryCondition<Type>& bc
)
:
    outflowBoundaryCondition<Type,staggered>(bc.mshField(), bc.mshBoundary())
{}

template<class Type>
stagOutflowBoundaryCondition<Type>::stagOutflowBoundaryCondition
(
    const meshField<Type,staggered>& field,
    const stagOutflowBoundaryCondition<Type>& bc
)
:
    outflowBoundaryCondition<Type,staggered>(field, bc.mshBoundary())
{}

template<class Type>
void stagOutflowBoundaryCondition<Type>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    // Eliminate only for non-shifted boundaries

    if (faceNumber(this->offset())/2 != d)
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
}

}

}

}

