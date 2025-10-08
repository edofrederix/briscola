#include "boundaryConditions.H"
#include "slipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(slip,slip,scalar,staggered)
addBoundaryConditionType(slip,scalar,staggered)

}

}

}
