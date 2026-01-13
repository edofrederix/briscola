#include "boundaryConditions.H"
#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(Neumann,Neumann,scalar,colocated)
makeBoundaryConditionType(Neumann,Neumann,vector,colocated)
makeBoundaryConditionType(Neumann,Neumann,tensor,colocated)

}

}

}
