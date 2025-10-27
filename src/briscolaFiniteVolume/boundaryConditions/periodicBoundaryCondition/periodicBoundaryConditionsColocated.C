#include "boundaryConditions.H"
#include "periodicBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the periodic boundary condition type for all possible data types

makeBoundaryConditionType(periodic,periodic,label,colocated)
makeBoundaryConditionType(periodic,periodic,scalar,colocated)
makeBoundaryConditionType(periodic,periodic,faceScalar,colocated)
makeBoundaryConditionType(periodic,periodic,edgeScalar,colocated)
makeBoundaryConditionType(periodic,periodic,vertexScalar,colocated)
makeBoundaryConditionType(periodic,periodic,vector,colocated)
makeBoundaryConditionType(periodic,periodic,faceVector,colocated)
makeBoundaryConditionType(periodic,periodic,edgeVector,colocated)
makeBoundaryConditionType(periodic,periodic,vertexVector,colocated)
makeBoundaryConditionType(periodic,periodic,tensor,colocated)
makeBoundaryConditionType(periodic,periodic,sphericalTensor,colocated)
makeBoundaryConditionType(periodic,periodic,symmTensor,colocated)
makeBoundaryConditionType(periodic,periodic,diagTensor,colocated)
makeBoundaryConditionType(periodic,periodic,stencil,colocated)
makeBoundaryConditionType(periodic,periodic,diagStencil,colocated)

}

}

}
