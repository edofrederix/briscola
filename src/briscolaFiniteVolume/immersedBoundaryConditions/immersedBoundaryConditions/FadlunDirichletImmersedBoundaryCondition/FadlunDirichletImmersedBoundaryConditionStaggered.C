#include "immersedBoundaryConditions.H"
#include "FadlunDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(FadlunDirichlet,scalar,staggered);
makeImmersedBoundaryConditionType(FadlunDirichlet,vector,staggered);

}

}

}
