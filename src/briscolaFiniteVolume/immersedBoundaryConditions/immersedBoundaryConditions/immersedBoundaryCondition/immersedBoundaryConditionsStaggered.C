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

makeImmersedBoundaryConditionBase(label,staggered);
makeImmersedBoundaryConditionBase(scalar,staggered);
makeImmersedBoundaryConditionBase(faceScalar,staggered);
makeImmersedBoundaryConditionBase(vertexScalar,staggered);
makeImmersedBoundaryConditionBase(vector,staggered);
makeImmersedBoundaryConditionBase(faceVector,staggered);
makeImmersedBoundaryConditionBase(vertexVector,staggered);
makeImmersedBoundaryConditionBase(tensor,staggered);
makeImmersedBoundaryConditionBase(symmTensor,staggered);
makeImmersedBoundaryConditionBase(stencil,staggered);
makeImmersedBoundaryConditionBase(diagStencil,staggered);

}

}

}
