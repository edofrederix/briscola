#include "boundaryConditions.H"
#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(Neumann,Neumann,staggered)
addBoundaryConditionTypes(Neumann,staggered)

}

}

}
