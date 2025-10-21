#include "DirichletBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void DirichletBoundaryConditionBase<Type,MeshType>::init
(
    const List<Type>& values
)
{
    boundaryValues_.clear();
    boundaryValues_.resize(this->fvMsh_.size()*MeshType::numberOfDirections);

    label item = 0;
    forAll(this->fvMsh_, l)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            boundaryValues_.set
            (
                item++,
                new block<Type>(this->N(l,d), values[d])
            );
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(List<Type>(this->dict().lookup("values")));
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& values
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(values);
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const FastPtrList<block<Type>>& boundaryValues
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

