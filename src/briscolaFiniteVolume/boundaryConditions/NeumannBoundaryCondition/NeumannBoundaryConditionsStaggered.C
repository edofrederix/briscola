#include "boundaryConditions.H"
#include "stagNeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Base
makeBoundaryConditionTypes(Neumann,Neumann,staggered)

// Derived
makeBoundaryConditionMeshedTypes(stagNeumann,Neumann,staggered)
addBoundaryConditionMeshedTypes(stagNeumann,staggered)

}

}

}
