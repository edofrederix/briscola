#include "immersedBoundaryConditions.H"
#include "penalizationDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(penalizationDirichlet,scalar,colocated);
makeImmersedBoundaryConditionType(penalizationDirichlet,vector,colocated);
makeImmersedBoundaryConditionType(penalizationDirichlet,tensor,colocated);

}

}

}
