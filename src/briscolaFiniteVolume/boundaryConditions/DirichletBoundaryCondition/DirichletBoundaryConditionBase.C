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
    const physicalType& value
)
{
    boundaryValues_.clear();
    boundaryValues_.resize(this->fvMsh_.size()*MeshType::numberOfDirections);

    tensor base(eye);

    if (this->fvMsh_.msh().template castable<rectilinearMesh>())
    {
        base = this->fvMsh_.msh().template cast<rectilinearMesh>().base();
    }

    label item = 0;
    forAll(this->fvMsh_, l)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            boundaryValues_.set
            (
                item++,
                new block<Type>
                (
                    this->N(l,d),
                    MeshType::project(value, d, base)
                )
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
    init(this->dict().template lookup<physicalType>("value"));
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const zero&
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(pTraits<physicalType>::zero);
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const physicalType& value
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(value);
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

