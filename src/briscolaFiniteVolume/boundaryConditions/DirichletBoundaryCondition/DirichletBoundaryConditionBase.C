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
    boundaryValues_.resize(MeshType::numberOfDirections);

    tensor base(eye);

    if (this->fvMsh_.msh().template castable<rectilinearMesh>())
    {
        base = this->fvMsh_.msh().template cast<rectilinearMesh>().base();
    }

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        boundaryValues_.set
        (
            d,
            new block<Type>
            (
                this->N(d),
                MeshType::project(value, d, base)
            )
        );
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(level, b)
{
    init(this->dict().template lookup<physicalType>("value"));
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b,
    const zero&
)
:
    boundaryCondition<Type,MeshType>(level, b)
{
    init(pTraits<physicalType>::zero);
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b,
    const physicalType& value
)
:
    boundaryCondition<Type,MeshType>(level, b)
{
    init(value);
}

template<class Type, class MeshType>
DirichletBoundaryConditionBase<Type,MeshType>::DirichletBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b,
    const FastPtrList<block<Type>>& boundaryValues
)
:
    boundaryCondition<Type,MeshType>(level, b),
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
    const DirichletBoundaryConditionBase<Type,MeshType>& bc,
    const meshLevel<Type,MeshType>& level
)
:
    boundaryCondition<Type,MeshType>(bc, level),
    boundaryValues_(bc.boundaryValues_)
{}

}

}

}
