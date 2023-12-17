#include "boundaryConditions.H"
#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionTypes(Dirichlet,Dirichlet,colocated)
addBoundaryConditionTypes(Dirichlet,colocated)

}

}

}
