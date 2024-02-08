#include "emptyBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    NeumannBoundaryCondition<Type,MeshType>
    (
        mshField,
        b,
        List<Type>(MeshType::numberOfDirections, Zero)
    )
{}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    NeumannBoundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary())
{}

template<class Type, class MeshType>
emptyBoundaryCondition<Type,MeshType>::emptyBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const emptyBoundaryCondition<Type,MeshType>& bc
)
:
    NeumannBoundaryCondition<Type,MeshType>(field, bc.mshBoundary())
{}

}

}

}

