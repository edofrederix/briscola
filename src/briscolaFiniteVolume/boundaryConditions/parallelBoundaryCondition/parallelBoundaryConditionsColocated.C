#include "boundaryConditions.H"
#include "parallelBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the parallel boundary condition type for all possible data types

makeBoundaryConditionType(parallel,parallel,label,colocated)
makeBoundaryConditionType(parallel,parallel,scalar,colocated)
makeBoundaryConditionType(parallel,parallel,faceScalar,colocated)
makeBoundaryConditionType(parallel,parallel,edgeScalar,colocated)
makeBoundaryConditionType(parallel,parallel,vertexScalar,colocated)
makeBoundaryConditionType(parallel,parallel,vector,colocated)
makeBoundaryConditionType(parallel,parallel,faceVector,colocated)
makeBoundaryConditionType(parallel,parallel,edgeVector,colocated)
makeBoundaryConditionType(parallel,parallel,vertexVector,colocated)
makeBoundaryConditionType(parallel,parallel,tensor,colocated)
makeBoundaryConditionType(parallel,parallel,sphericalTensor,colocated)
makeBoundaryConditionType(parallel,parallel,symmTensor,colocated)
makeBoundaryConditionType(parallel,parallel,diagTensor,colocated)
makeBoundaryConditionType(parallel,parallel,stencil,colocated)
makeBoundaryConditionType(parallel,parallel,diagStencil,colocated)

}

}

}
