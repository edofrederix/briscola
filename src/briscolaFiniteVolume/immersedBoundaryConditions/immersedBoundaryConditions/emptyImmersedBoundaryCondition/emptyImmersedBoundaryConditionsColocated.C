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

makeImmersedBoundaryConditionType(empty,label,colocated);
makeImmersedBoundaryConditionType(empty,scalar,colocated);
makeImmersedBoundaryConditionType(empty,faceScalar,colocated);
makeImmersedBoundaryConditionType(empty,edgeScalar,colocated);
makeImmersedBoundaryConditionType(empty,vertexScalar,colocated);
makeImmersedBoundaryConditionType(empty,vector,colocated);
makeImmersedBoundaryConditionType(empty,faceVector,colocated);
makeImmersedBoundaryConditionType(empty,edgeVector,colocated);
makeImmersedBoundaryConditionType(empty,vertexVector,colocated);
makeImmersedBoundaryConditionType(empty,tensor,colocated);
makeImmersedBoundaryConditionType(empty,sphericalTensor,colocated)
makeImmersedBoundaryConditionType(empty,symmTensor,colocated);
makeImmersedBoundaryConditionType(empty,diagTensor,colocated);
makeImmersedBoundaryConditionType(empty,stencil,colocated);
makeImmersedBoundaryConditionType(empty,diagStencil,colocated);

}

}

}
