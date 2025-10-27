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
    boundaryGradients_.resize(this->fvMsh_.size()*MeshType::numberOfDirections);

    tensor base(eye);

    if (this->fvMsh_.msh().template castable<rectilinearMesh>())
    {
        base = this->fvMsh_.msh().template cast<rectilinearMesh>().base();
    }

    label item = 0;
    forAll(this->fvMsh_, l)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            boundaryGradients_.set
            (
                item++,
                new block<Type>
                (
                    this->N(l,d),
                    MeshType::project(gradient, d, base)
                )
            );
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(this->dict().template lookup<physicalType>("gradient"));
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
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
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const physicalType& gradient
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(gradient);
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const FastPtrList<block<Type>>& boundaryGradients
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

