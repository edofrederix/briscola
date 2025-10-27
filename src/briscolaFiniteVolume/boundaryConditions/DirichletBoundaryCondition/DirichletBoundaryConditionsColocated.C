#include "boundaryConditions.H"
#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(Dirichlet,Dirichlet,scalar,colocated)
makeBoundaryConditionType(Dirichlet,Dirichlet,vector,colocated)
makeBoundaryConditionType(Dirichlet,Dirichlet,tensor,colocated)

}

}

}
