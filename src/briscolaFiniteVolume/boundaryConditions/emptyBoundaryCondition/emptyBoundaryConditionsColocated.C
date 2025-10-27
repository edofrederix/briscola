#include "boundaryConditions.H"
#include "emptyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the empty boundary condition type for all possible data types

makeBoundaryConditionType(empty,empty,label,colocated)
makeBoundaryConditionType(empty,empty,scalar,colocated)
makeBoundaryConditionType(empty,empty,faceScalar,colocated)
makeBoundaryConditionType(empty,empty,edgeScalar,colocated)
makeBoundaryConditionType(empty,empty,vertexScalar,colocated)
makeBoundaryConditionType(empty,empty,vector,colocated)
makeBoundaryConditionType(empty,empty,faceVector,colocated)
makeBoundaryConditionType(empty,empty,edgeVector,colocated)
makeBoundaryConditionType(empty,empty,vertexVector,colocated)
makeBoundaryConditionType(empty,empty,tensor,colocated)
makeBoundaryConditionType(empty,empty,sphericalTensor,colocated)
makeBoundaryConditionType(empty,empty,symmTensor,colocated)
makeBoundaryConditionType(empty,empty,diagTensor,colocated)
makeBoundaryConditionType(empty,empty,stencil,colocated)
makeBoundaryConditionType(empty,empty,diagStencil,colocated)

}

}

}
