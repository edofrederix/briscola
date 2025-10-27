#include "emptyBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
emptyBoundaryConditionBase<Type,MeshType>::emptyBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{}

template<class Type, class MeshType>
emptyBoundaryConditionBase<Type,MeshType>::emptyBoundaryConditionBase
(
    const emptyBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc)
{}

template<class Type, class MeshType>
emptyBoundaryConditionBase<Type,MeshType>::emptyBoundaryConditionBase
(
    const meshField<Type,MeshType>& field,
    const emptyBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc)
{}

}

}

}

