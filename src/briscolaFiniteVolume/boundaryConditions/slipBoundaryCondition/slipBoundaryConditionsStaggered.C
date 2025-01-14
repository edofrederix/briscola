#include "boundaryConditions.H"
#include "stagSlipBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Derived
makeBoundaryConditionMeshedType(stagSlip,slip,scalar,staggered)
addBoundaryConditionMeshedType(stagSlip,scalar,staggered)

}

}

}
