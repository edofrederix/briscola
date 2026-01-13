#include "boundaryConditions.H"
#include "noSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(noSlip,noSlip,scalar,colocated)
makeBoundaryConditionType(noSlip,noSlip,vector,colocated)
makeBoundaryConditionType(noSlip,noSlip,tensor,colocated)

}

}

}
