#include "boundaryConditions.H"
#include "dummyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the dummy boundary condition type for all possible data types

makeBoundaryConditionType(dummy,dummy,label,colocated)
makeBoundaryConditionType(dummy,dummy,scalar,colocated)
makeBoundaryConditionType(dummy,dummy,faceScalar,colocated)
makeBoundaryConditionType(dummy,dummy,edgeScalar,colocated)
makeBoundaryConditionType(dummy,dummy,vertexScalar,colocated)
makeBoundaryConditionType(dummy,dummy,vector,colocated)
makeBoundaryConditionType(dummy,dummy,faceVector,colocated)
makeBoundaryConditionType(dummy,dummy,edgeVector,colocated)
makeBoundaryConditionType(dummy,dummy,vertexVector,colocated)
makeBoundaryConditionType(dummy,dummy,tensor,colocated)
makeBoundaryConditionType(dummy,dummy,sphericalTensor,colocated)
makeBoundaryConditionType(dummy,dummy,symmTensor,colocated)
makeBoundaryConditionType(dummy,dummy,diagTensor,colocated)
makeBoundaryConditionType(dummy,dummy,stencil,colocated)
makeBoundaryConditionType(dummy,dummy,diagStencil,colocated)

}

}

}
