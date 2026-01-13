#include "immersedBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// We need the empty immersed boundary condition type for all possible data
// types because it is used by default

makeImmersedBoundaryConditionType(empty,label,staggered);
makeImmersedBoundaryConditionType(empty,scalar,staggered);
makeImmersedBoundaryConditionType(empty,faceScalar,staggered);
makeImmersedBoundaryConditionType(empty,edgeScalar,staggered);
makeImmersedBoundaryConditionType(empty,vertexScalar,staggered);
makeImmersedBoundaryConditionType(empty,vector,staggered);
makeImmersedBoundaryConditionType(empty,faceVector,staggered);
makeImmersedBoundaryConditionType(empty,edgeVector,staggered);
makeImmersedBoundaryConditionType(empty,vertexVector,staggered);
makeImmersedBoundaryConditionType(empty,tensor,staggered);
makeImmersedBoundaryConditionType(empty,sphericalTensor,staggered)
makeImmersedBoundaryConditionType(empty,symmTensor,staggered);
makeImmersedBoundaryConditionType(empty,diagTensor,staggered);
makeImmersedBoundaryConditionType(empty,stencil,staggered);
makeImmersedBoundaryConditionType(empty,diagStencil,staggered);

}

}

}
