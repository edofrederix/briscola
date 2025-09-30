#include "boundaryConditions.H"
#include "mappedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(mapped,mapped,scalar,staggered)
makeBoundaryConditionType(mapped,mapped,vector,staggered)
makeBoundaryConditionType(mapped,mapped,tensor,staggered)
makeBoundaryConditionType(mapped,mapped,sphericalTensor,staggered)
makeBoundaryConditionType(mapped,mapped,symmTensor,staggered)
makeBoundaryConditionType(mapped,mapped,diagTensor,staggered)
makeBoundaryConditionType(mapped,mapped,faceScalar,staggered)
makeBoundaryConditionType(mapped,mapped,faceVector,staggered)

addBoundaryConditionType(mapped,scalar,staggered)
addBoundaryConditionType(mapped,vector,staggered)
addBoundaryConditionType(mapped,tensor,staggered)
addBoundaryConditionType(mapped,sphericalTensor,staggered)
addBoundaryConditionType(mapped,symmTensor,staggered)
addBoundaryConditionType(mapped,diagTensor,staggered)
addBoundaryConditionType(mapped,faceScalar,staggered)
addBoundaryConditionType(mapped,faceVector,staggered)

}

}

}
