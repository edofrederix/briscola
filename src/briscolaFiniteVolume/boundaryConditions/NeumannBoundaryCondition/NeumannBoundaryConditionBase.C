#include "NeumannBoundaryConditionBase.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryGradients_(this->dict().lookup("gradients"))
{}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& boundaryGradients
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryGradients_(boundaryGradients)
{}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const NeumannBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc),
    boundaryGradients_(bc.boundaryGradients_)
{}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& field,
    const NeumannBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc),
    boundaryGradients_(bc.boundaryGradients_)
{}

}

}

}

