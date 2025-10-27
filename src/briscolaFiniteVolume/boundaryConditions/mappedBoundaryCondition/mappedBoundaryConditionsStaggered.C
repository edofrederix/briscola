#include "boundaryConditions.H"
#include "mappedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(mapped,mapped,scalar,staggered)
makeBoundaryConditionType(mapped,mapped,vector,staggered)

}

}

}
