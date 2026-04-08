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
    patchBoundaryCondition<Type,MeshType>(level, b)
{}

template<class Type, class MeshType>
emptyBoundaryConditionBase<Type,MeshType>::emptyBoundaryConditionBase
(
    const emptyBoundaryConditionBase<Type,MeshType>& bc
)
:
    patchBoundaryCondition<Type,MeshType>(bc)
{}

template<class Type, class MeshType>
emptyBoundaryConditionBase<Type,MeshType>::emptyBoundaryConditionBase
(
    const emptyBoundaryConditionBase<Type,MeshType>& bc,
    const meshLevel<Type,MeshType>& level
)
:
    patchBoundaryCondition<Type,MeshType>(bc, level)
{}

}

}

}

