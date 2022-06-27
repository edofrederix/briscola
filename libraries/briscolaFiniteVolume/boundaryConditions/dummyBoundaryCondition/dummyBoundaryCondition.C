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

makeBoundaryConditionType(dummy,label,colocated)
makeBoundaryConditionType(dummy,label,staggered)

makeBoundaryConditionType(dummy,scalar,colocated)
makeBoundaryConditionType(dummy,scalar,staggered)

makeBoundaryConditionType(dummy,faceScalar,colocated)
makeBoundaryConditionType(dummy,faceScalar,staggered)

makeBoundaryConditionType(dummy,vector,colocated)
makeBoundaryConditionType(dummy,vector,staggered)

makeBoundaryConditionType(dummy,faceVector,colocated)
makeBoundaryConditionType(dummy,faceVector,staggered)

makeBoundaryConditionType(dummy,tensor,colocated)
makeBoundaryConditionType(dummy,tensor,staggered)

makeBoundaryConditionType(dummy,sphericalTensor,colocated)
makeBoundaryConditionType(dummy,sphericalTensor,staggered)

makeBoundaryConditionType(dummy,symmTensor,colocated)
makeBoundaryConditionType(dummy,symmTensor,staggered)

makeBoundaryConditionType(dummy,diagTensor,colocated)
makeBoundaryConditionType(dummy,diagTensor,staggered)

makeBoundaryConditionType(dummy,stencil,colocated)
makeBoundaryConditionType(dummy,stencil,staggered)

makeBoundaryConditionType(dummy,diagStencil,colocated)
makeBoundaryConditionType(dummy,diagStencil,staggered)

}

}

}

