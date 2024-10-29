#include "boundaryConditions.H"
#include "stagOutflowBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Base
makeBoundaryConditionType(outflow,outflow,scalar,staggered)

// Derived
makeBoundaryConditionMeshedType(stagOutflow,outflow,scalar,staggered)
addBoundaryConditionMeshedType(stagOutflow,scalar,staggered)

}

}

}
