#include "boundaryConditions.H"
#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(Neumann,Neumann,scalar,staggered)
makeBoundaryConditionType(Neumann,Neumann,vector,staggered)

}

}

}
