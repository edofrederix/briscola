#include "boundaryConditions.H"
#include "periodicBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the periodic boundary condition type for all possible data types

makeBoundaryConditionType(periodic,periodic,label,staggered)
makeBoundaryConditionType(periodic,periodic,scalar,staggered)
makeBoundaryConditionType(periodic,periodic,faceScalar,staggered)
makeBoundaryConditionType(periodic,periodic,edgeScalar,staggered)
makeBoundaryConditionType(periodic,periodic,vertexScalar,staggered)
makeBoundaryConditionType(periodic,periodic,vector,staggered)
makeBoundaryConditionType(periodic,periodic,faceVector,staggered)
makeBoundaryConditionType(periodic,periodic,edgeVector,staggered)
makeBoundaryConditionType(periodic,periodic,vertexVector,staggered)
makeBoundaryConditionType(periodic,periodic,tensor,staggered)
makeBoundaryConditionType(periodic,periodic,sphericalTensor,staggered)
makeBoundaryConditionType(periodic,periodic,symmTensor,staggered)
makeBoundaryConditionType(periodic,periodic,diagTensor,staggered)
makeBoundaryConditionType(periodic,periodic,stencil,staggered)
makeBoundaryConditionType(periodic,periodic,diagStencil,staggered)

}

}

}
