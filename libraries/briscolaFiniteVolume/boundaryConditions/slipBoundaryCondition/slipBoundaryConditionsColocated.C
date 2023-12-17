#include "boundaryConditions.H"
#include "slipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(slip,slip,vector,colocated)
addBoundaryConditionType(slip,vector,colocated)

}

}

}
