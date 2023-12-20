#include "noSlipBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const boundary& b
)
:
    DirichletBoundaryCondition<Type,MeshType>
    (
        mshField,
        b,
        List<Type>(MeshType::numberOfDirections, Zero)
    )
{}

template<class Type, class MeshType>
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const noSlipBoundaryCondition<Type,MeshType>& bc
)
:
    DirichletBoundaryCondition<Type,MeshType>(bc.mshField(), bc.mshBoundary())
{}

template<class Type, class MeshType>
noSlipBoundaryCondition<Type,MeshType>::noSlipBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const noSlipBoundaryCondition<Type,MeshType>& bc
)
:
    DirichletBoundaryCondition<Type,MeshType>(field, bc.mshBoundary())
{}

}

}

}

