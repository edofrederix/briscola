#include "immersedBoundaryConditions.H"
#include "MittalNeumannImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeImmersedBoundaryConditionType(MittalNeumann,scalar,staggered);
makeImmersedBoundaryConditionType(MittalNeumann,vector,staggered);

}

}

}
