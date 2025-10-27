#include "boundaryConditions.H"
#include "mappedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(mapped,mapped,scalar,colocated)
makeBoundaryConditionType(mapped,mapped,vector,colocated)
makeBoundaryConditionType(mapped,mapped,tensor,colocated)

}

}

}
