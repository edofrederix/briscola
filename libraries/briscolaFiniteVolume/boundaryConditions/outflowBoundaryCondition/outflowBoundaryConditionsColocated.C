#include "boundaryConditions.H"
#include "outflowBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(outflow,outflow,vector,colocated)
addBoundaryConditionType(outflow,vector,colocated)

}

}

}
