#include "immersedBoundaryConditions.H"
#include "penalizationDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(penalizationDirichlet,scalar,staggered);
makeImmersedBoundaryConditionType(penalizationDirichlet,vector,staggered);

}

}

}
