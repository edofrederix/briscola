#include "boundaryConditions.H"
#include "noSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(noSlip,noSlip,vector,colocated)
addBoundaryConditionType(noSlip,vector,colocated)

}

}

}
