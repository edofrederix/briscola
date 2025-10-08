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
    const List<Type>& gradients
)
{
    boundaryGradients_.clear();
    boundaryGradients_.resize(this->fvMsh_.size()*MeshType::numberOfDirections);

    label c = 0;
    forAll(this->fvMsh_, l)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            boundaryGradients_.set
            (
                c++,
                new block<Type>(this->N(l,d), gradients[d])
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
    init(List<Type>(this->dict().lookup("gradients")));
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const List<Type>& gradients
)
:
    boundaryCondition<Type,MeshType>(mshField, b)
{
    init(gradients);
}

template<class Type, class MeshType>
NeumannBoundaryConditionBase<Type,MeshType>::NeumannBoundaryConditionBase
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b,
    const PtrList<block<Type>>& boundaryGradients
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

