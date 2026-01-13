#include "immersedBoundaryConditions.H"
#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(VremanDirichlet,scalar,staggered);
makeImmersedBoundaryConditionType(VremanDirichlet,vector,staggered);

}

}

}
