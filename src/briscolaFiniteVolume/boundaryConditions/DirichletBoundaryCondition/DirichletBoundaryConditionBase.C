#include "DirichletBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryValues_(this->dict().lookup("values"))
{}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& boundaryValues
)
:
    boundaryCondition<Type,MeshType>(mshField, b),
    boundaryValues_(boundaryValues)
{}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const DirichletBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc),
    boundaryValues_(bc.boundaryValues_)
{}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& field,
    const DirichletBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc),
    boundaryValues_(bc.boundaryValues_)
{}

}

}

}

