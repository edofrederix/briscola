#include "immersedBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

#include "emptyImmersedBoundaryCondition.H"
#include "penalizationDirichletImmersedBoundaryCondition.H"
#include "FadlunDirichletImmersedBoundaryCondition.H"
#include "MittalDirichletImmersedBoundaryCondition.H"
#include "MittalNeumannImmersedBoundaryCondition.H"
#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionBase(label,colocated);
makeImmersedBoundaryConditionBase(scalar,colocated);
makeImmersedBoundaryConditionBase(faceScalar,colocated);
makeImmersedBoundaryConditionBase(edgeScalar,colocated);
makeImmersedBoundaryConditionBase(vertexScalar,colocated);
makeImmersedBoundaryConditionBase(vector,colocated);
makeImmersedBoundaryConditionBase(faceVector,colocated);
makeImmersedBoundaryConditionBase(edgeVector,colocated);
makeImmersedBoundaryConditionBase(vertexVector,colocated);
makeImmersedBoundaryConditionBase(tensor,colocated);
makeImmersedBoundaryConditionBase(sphericalTensor,colocated)
makeImmersedBoundaryConditionBase(symmTensor,colocated);
makeImmersedBoundaryConditionBase(diagTensor,colocated);
makeImmersedBoundaryConditionBase(stencil,colocated);
makeImmersedBoundaryConditionBase(diagStencil,colocated);

}

}

}
