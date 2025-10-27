#include "boundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionBase(label,colocated)
makeBoundaryConditionBase(scalar,colocated)
makeBoundaryConditionBase(faceScalar,colocated)
makeBoundaryConditionBase(edgeScalar,colocated)
makeBoundaryConditionBase(vertexScalar,colocated)
makeBoundaryConditionBase(vector,colocated)
makeBoundaryConditionBase(faceVector,colocated)
makeBoundaryConditionBase(edgeVector,colocated)
makeBoundaryConditionBase(vertexVector,colocated)
makeBoundaryConditionBase(tensor,colocated)
makeBoundaryConditionBase(sphericalTensor,colocated)
makeBoundaryConditionBase(symmTensor,colocated)
makeBoundaryConditionBase(diagTensor,colocated)
makeBoundaryConditionBase(stencil,colocated)
makeBoundaryConditionBase(diagStencil,colocated)

}

}

}
