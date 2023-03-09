#include "dummyBoundaryCondition.H"

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
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const partPatch& patch
)
:
    boundaryCondition<Type,MeshType>(mshField, patch)
{}

template<class Type, class MeshType>
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const dummyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(bc.mshField(), bc.patch())
{}

template<class Type, class MeshType>
dummyBoundaryCondition<Type,MeshType>::dummyBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const dummyBoundaryCondition<Type,MeshType>& bc
)
:
    boundaryCondition<Type,MeshType>(field, bc.patch())
{}

template<class Type, class MeshType>
void dummyBoundaryCondition<Type,MeshType>::initEvaluate(const label)
{}

template<class Type, class MeshType>
void dummyBoundaryCondition<Type,MeshType>::evaluate
(
    const label,
    const bool
)
{}

makeBoundaryConditionTypes(dummy)


}

}

}

