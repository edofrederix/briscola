#include "immersedBoundaryConditions.H"
#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(MittalDirichlet,scalar,colocated);
makeImmersedBoundaryConditionType(MittalDirichlet,vector,colocated);
makeImmersedBoundaryConditionType(MittalDirichlet,tensor,colocated);

}

}

}
