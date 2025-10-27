#include "immersedBoundaryConditions.H"
#include "MittalNeumannImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(MittalNeumann,scalar,colocated);
makeImmersedBoundaryConditionType(MittalNeumann,vector,colocated);
makeImmersedBoundaryConditionType(MittalNeumann,tensor,colocated);

}

}

}
