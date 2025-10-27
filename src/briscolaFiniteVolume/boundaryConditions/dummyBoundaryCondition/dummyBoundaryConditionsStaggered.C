#include "boundaryConditions.H"
#include "dummyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the dummy boundary condition type for all possible data types

makeBoundaryConditionType(dummy,dummy,label,staggered)
makeBoundaryConditionType(dummy,dummy,scalar,staggered)
makeBoundaryConditionType(dummy,dummy,faceScalar,staggered)
makeBoundaryConditionType(dummy,dummy,edgeScalar,staggered)
makeBoundaryConditionType(dummy,dummy,vertexScalar,staggered)
makeBoundaryConditionType(dummy,dummy,vector,staggered)
makeBoundaryConditionType(dummy,dummy,faceVector,staggered)
makeBoundaryConditionType(dummy,dummy,edgeVector,staggered)
makeBoundaryConditionType(dummy,dummy,vertexVector,staggered)
makeBoundaryConditionType(dummy,dummy,tensor,staggered)
makeBoundaryConditionType(dummy,dummy,sphericalTensor,staggered)
makeBoundaryConditionType(dummy,dummy,symmTensor,staggered)
makeBoundaryConditionType(dummy,dummy,diagTensor,staggered)
makeBoundaryConditionType(dummy,dummy,stencil,staggered)
makeBoundaryConditionType(dummy,dummy,diagStencil,staggered)

}

}

}
