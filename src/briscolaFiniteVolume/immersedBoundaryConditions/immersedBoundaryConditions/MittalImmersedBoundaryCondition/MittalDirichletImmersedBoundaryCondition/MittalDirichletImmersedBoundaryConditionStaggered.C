#include "immersedBoundaryConditions.H"
#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(MittalDirichlet,scalar,staggered);
makeImmersedBoundaryConditionType(MittalDirichlet,vector,staggered);

}

}

}
