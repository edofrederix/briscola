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
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(level, b)
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
    const emptyBoundaryConditionBase<Type,MeshType>& bc,
    const meshLevel<Type,MeshType>& level
)
:
    boundaryCondition<Type,MeshType>(bc, level)
{}

}

}

}

