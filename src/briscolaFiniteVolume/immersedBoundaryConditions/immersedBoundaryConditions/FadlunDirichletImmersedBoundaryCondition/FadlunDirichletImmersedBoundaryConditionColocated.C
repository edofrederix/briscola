#include "immersedBoundaryConditions.H"
#include "FadlunDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(FadlunDirichlet,scalar,colocated);
makeImmersedBoundaryConditionType(FadlunDirichlet,vector,colocated);
makeImmersedBoundaryConditionType(FadlunDirichlet,tensor,colocated);

}

}

}
