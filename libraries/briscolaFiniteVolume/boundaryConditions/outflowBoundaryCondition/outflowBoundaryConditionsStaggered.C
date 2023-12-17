#include "boundaryConditions.H"
#include "stagOutflowBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Base
makeBoundaryConditionTypes(outflow,outflow,staggered)

// Derived
makeBoundaryConditionMeshedTypes(stagOutflow,outflow,staggered)
addBoundaryConditionMeshedTypes(stagOutflow,staggered)

}

}

}
