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
makeBoundaryConditionBase(vertexScalar,colocated)
makeBoundaryConditionBase(vector,colocated)
makeBoundaryConditionBase(faceVector,colocated)
makeBoundaryConditionBase(vertexVector,colocated)
makeBoundaryConditionBase(tensor,colocated)
makeBoundaryConditionBase(symmTensor,colocated)
makeBoundaryConditionBase(stencil,colocated)
makeBoundaryConditionBase(diagStencil,colocated)

}

}

}
