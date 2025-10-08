#include "boundaryConditions.H"
#include "noSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(noSlip,noSlip,scalar,staggered)
addBoundaryConditionType(noSlip,scalar,staggered)

}

}

}
