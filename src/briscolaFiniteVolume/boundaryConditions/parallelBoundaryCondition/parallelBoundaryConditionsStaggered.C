#include "boundaryConditions.H"
#include "parallelBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the parallel boundary condition type for all possible data types

makeBoundaryConditionType(parallel,parallel,label,staggered)
makeBoundaryConditionType(parallel,parallel,scalar,staggered)
makeBoundaryConditionType(parallel,parallel,faceScalar,staggered)
makeBoundaryConditionType(parallel,parallel,edgeScalar,staggered)
makeBoundaryConditionType(parallel,parallel,vertexScalar,staggered)
makeBoundaryConditionType(parallel,parallel,vector,staggered)
makeBoundaryConditionType(parallel,parallel,faceVector,staggered)
makeBoundaryConditionType(parallel,parallel,edgeVector,staggered)
makeBoundaryConditionType(parallel,parallel,vertexVector,staggered)
makeBoundaryConditionType(parallel,parallel,tensor,staggered)
makeBoundaryConditionType(parallel,parallel,sphericalTensor,staggered)
makeBoundaryConditionType(parallel,parallel,symmTensor,staggered)
makeBoundaryConditionType(parallel,parallel,diagTensor,staggered)
makeBoundaryConditionType(parallel,parallel,stencil,staggered)
makeBoundaryConditionType(parallel,parallel,diagStencil,staggered)

}

}

}
