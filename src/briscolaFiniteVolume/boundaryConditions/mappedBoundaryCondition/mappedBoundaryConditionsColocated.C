#include "boundaryConditions.H"
#include "mappedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Types supported by data exchange class

makeBoundaryConditionType(mapped,mapped,scalar,colocated)
makeBoundaryConditionType(mapped,mapped,vector,colocated)
makeBoundaryConditionType(mapped,mapped,tensor,colocated)
makeBoundaryConditionType(mapped,mapped,sphericalTensor,colocated)
makeBoundaryConditionType(mapped,mapped,symmTensor,colocated)
makeBoundaryConditionType(mapped,mapped,diagTensor,colocated)
makeBoundaryConditionType(mapped,mapped,faceScalar,colocated)
makeBoundaryConditionType(mapped,mapped,faceVector,colocated)

addBoundaryConditionType(mapped,scalar,colocated)
addBoundaryConditionType(mapped,vector,colocated)
addBoundaryConditionType(mapped,tensor,colocated)
addBoundaryConditionType(mapped,sphericalTensor,colocated)
addBoundaryConditionType(mapped,symmTensor,colocated)
addBoundaryConditionType(mapped,diagTensor,colocated)
addBoundaryConditionType(mapped,faceScalar,colocated)
addBoundaryConditionType(mapped,faceVector,colocated)

}

}

}
