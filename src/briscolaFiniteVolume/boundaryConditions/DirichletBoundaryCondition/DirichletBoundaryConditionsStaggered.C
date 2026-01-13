#include "boundaryConditions.H"
#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionType(Dirichlet,Dirichlet,scalar,staggered)
makeBoundaryConditionType(Dirichlet,Dirichlet,vector,staggered)

}

}

}
