#include "boundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

#include "DirichletBoundaryCondition.H"
#include "NeumannBoundaryCondition.H"

#include "parallelBoundaryCondition.H"
#include "periodicBoundaryCondition.H"

#include "emptyBoundaryCondition.H"
#include "dummyBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Boundary condition types

makeBoundaryConditionTypes(Dirichlet)
makeBoundaryConditionTypes(Neumann)

makeBoundaryConditionTypes(parallel)
makeBoundaryConditionTypes(periodic)

makeBoundaryConditionTypes(empty)
makeBoundaryConditionTypes(dummy)

// Base types

makeBoundaryCondition(label,colocated)
makeBoundaryCondition(label,staggered)

makeBoundaryCondition(scalar,colocated)
makeBoundaryCondition(scalar,staggered)

makeBoundaryCondition(faceScalar,colocated)
makeBoundaryCondition(faceScalar,staggered)

makeBoundaryCondition(edgeScalar,colocated)
makeBoundaryCondition(edgeScalar,staggered)

makeBoundaryCondition(vertexScalar,colocated)
makeBoundaryCondition(vertexScalar,staggered)

makeBoundaryCondition(vector,colocated)
makeBoundaryCondition(vector,staggered)

makeBoundaryCondition(faceVector,colocated)
makeBoundaryCondition(faceVector,staggered)

makeBoundaryCondition(edgeVector,colocated)
makeBoundaryCondition(edgeVector,staggered)

makeBoundaryCondition(vertexVector,colocated)
makeBoundaryCondition(vertexVector,staggered)

makeBoundaryCondition(tensor,colocated)
makeBoundaryCondition(tensor,staggered)

makeBoundaryCondition(sphericalTensor,colocated)
makeBoundaryCondition(sphericalTensor,staggered)

makeBoundaryCondition(symmTensor,colocated)
makeBoundaryCondition(symmTensor,staggered)

makeBoundaryCondition(diagTensor,colocated)
makeBoundaryCondition(diagTensor,staggered)

makeBoundaryCondition(diagStencil,colocated)
makeBoundaryCondition(diagStencil,staggered)

makeBoundaryCondition(stencil,colocated)
makeBoundaryCondition(stencil,staggered)

}

}

}
