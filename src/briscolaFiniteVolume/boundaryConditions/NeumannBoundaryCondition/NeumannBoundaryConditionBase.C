#include "NeumannBoundaryConditionBase.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void NeumannBoundaryConditionBase<Type,MeshType>::init
(
    const physicalType& gradient
)
{
    boundaryGradients_.clear();
    boundaryGradients_.resize(MeshType::numberOfDirections);

    tensor base(eye);

    if (this->fvMsh_.msh().template castable<rectilinearMesh>())
    {
        base = this->fvMsh_.msh().template cast<rectilinearMesh>().base();
    }

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        boundaryGradients_.set
        (
            d,
            new block<Type>
            (
                this->N(d),
                MeshType::project(gradient, d, base)
            )
        );
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(level, b)
{
    init(this->dict().template lookup<physicalType>("gradient"));
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
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
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b,
    const physicalType& gradient
)
:
    boundaryCondition<Type,MeshType>(level, b)
{
    init(gradient);
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshLevel<Type,MeshType>& level,
    const boundary& b,
    const FastPtrList<block<Type>>& boundaryGradients
)
:
    boundaryCondition<Type,MeshType>(level, b),
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
    const meshLevel<Type,MeshType>& level,
    const NeumannBoundaryConditionBase<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(level, bc),
    boundaryGradients_(bc.boundaryGradients_)
{}

}

}

}

