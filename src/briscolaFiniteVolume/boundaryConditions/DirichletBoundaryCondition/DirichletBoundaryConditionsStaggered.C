#include "boundaryConditions.H"
#include "stagDirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Base
makeBoundaryConditionTypes(Dirichlet,Dirichlet,staggered)

// Derived
makeBoundaryConditionMeshedTypes(stagDirichlet,Dirichlet,staggered)
addBoundaryConditionMeshedTypes(stagDirichlet,staggered)

}

}

}
