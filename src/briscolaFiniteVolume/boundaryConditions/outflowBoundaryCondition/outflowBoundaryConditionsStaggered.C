#include "boundaryConditions.H"
#include "outflowBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(outflow,outflow,scalar,staggered)
addBoundaryConditionType(outflow,scalar,staggered)

}

}

}
