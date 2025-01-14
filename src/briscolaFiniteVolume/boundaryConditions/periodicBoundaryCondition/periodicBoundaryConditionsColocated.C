#include "boundaryConditions.H"
#include "periodicBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(periodic,periodic,colocated)
addBoundaryConditionTypes(periodic,colocated)

}

}

}
