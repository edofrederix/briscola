#include "boundaryConditions.H"
#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(Neumann,Neumann,colocated)
addBoundaryConditionTypes(Neumann,colocated)

}

}

}
