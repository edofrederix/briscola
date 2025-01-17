#include "boundaryConditions.H"
#include "parallelBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(parallel,parallel,staggered)
addBoundaryConditionTypes(parallel,staggered)

}

}

}
