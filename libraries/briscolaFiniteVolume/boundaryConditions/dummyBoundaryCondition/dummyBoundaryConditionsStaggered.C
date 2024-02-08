#include "boundaryConditions.H"
#include "dummyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(dummy,dummy,staggered)
addBoundaryConditionTypes(dummy,staggered)

}

}

}
