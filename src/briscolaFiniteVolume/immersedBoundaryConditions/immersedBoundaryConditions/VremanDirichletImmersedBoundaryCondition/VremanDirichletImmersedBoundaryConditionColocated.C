#include "immersedBoundaryConditions.H"
#include "VremanDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(VremanDirichlet,scalar,colocated);
makeImmersedBoundaryConditionType(VremanDirichlet,vector,colocated);
makeImmersedBoundaryConditionType(VremanDirichlet,tensor,colocated);

}

}

}
