#include "boundaryConditions.H"
#include "emptyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the empty boundary condition type for all possible data types

makeBoundaryConditionType(empty,empty,label,staggered)
makeBoundaryConditionType(empty,empty,scalar,staggered)
makeBoundaryConditionType(empty,empty,faceScalar,staggered)
makeBoundaryConditionType(empty,empty,edgeScalar,staggered)
makeBoundaryConditionType(empty,empty,vertexScalar,staggered)
makeBoundaryConditionType(empty,empty,vector,staggered)
makeBoundaryConditionType(empty,empty,faceVector,staggered)
makeBoundaryConditionType(empty,empty,edgeVector,staggered)
makeBoundaryConditionType(empty,empty,vertexVector,staggered)
makeBoundaryConditionType(empty,empty,tensor,staggered)
makeBoundaryConditionType(empty,empty,sphericalTensor,staggered)
makeBoundaryConditionType(empty,empty,symmTensor,staggered)
makeBoundaryConditionType(empty,empty,diagTensor,staggered)
makeBoundaryConditionType(empty,empty,stencil,staggered)
makeBoundaryConditionType(empty,empty,diagStencil,staggered)

}

}

}
