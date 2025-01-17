#include "boundaryConditions.H"
#include "emptyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(empty,empty,colocated)
addBoundaryConditionTypes(empty,colocated)

}

}

}
