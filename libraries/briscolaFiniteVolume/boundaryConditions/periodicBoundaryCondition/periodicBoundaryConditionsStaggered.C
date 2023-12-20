#include "boundaryConditions.H"
#include "periodicBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(periodic,periodic,staggered)
addBoundaryConditionTypes(periodic,staggered)

}

}

}
