#include "boundaryConditions.H"
#include "stagNoSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Derived
makeBoundaryConditionMeshedType(stagNoSlip,noSlip,scalar,staggered)
addBoundaryConditionMeshedType(stagNoSlip,scalar,staggered)

}

}

}
