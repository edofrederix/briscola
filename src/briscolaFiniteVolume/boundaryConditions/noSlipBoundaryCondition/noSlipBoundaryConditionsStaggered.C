#include "boundaryConditions.H"
#include "noSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(noSlip,noSlip,scalar,staggered)
makeBoundaryConditionType(noSlip,noSlip,vector,staggered)

}

}

}
