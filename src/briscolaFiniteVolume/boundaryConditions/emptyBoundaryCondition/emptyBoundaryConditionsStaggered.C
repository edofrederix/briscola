#include "boundaryConditions.H"
#include "stagEmptyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Base
makeBoundaryConditionTypes(empty,empty,staggered)

// Derived
makeBoundaryConditionMeshedTypes(stagEmpty,empty,staggered)
addBoundaryConditionMeshedTypes(stagEmpty,staggered)

}

}

}
