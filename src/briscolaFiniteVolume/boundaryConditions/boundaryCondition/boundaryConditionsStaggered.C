#include "boundaryConditions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeBoundaryConditionBase(label,staggered)
makeBoundaryConditionBase(scalar,staggered)
makeBoundaryConditionBase(faceScalar,staggered)
makeBoundaryConditionBase(edgeScalar,staggered)
makeBoundaryConditionBase(vertexScalar,staggered)
makeBoundaryConditionBase(vector,staggered)
makeBoundaryConditionBase(faceVector,staggered)
makeBoundaryConditionBase(edgeVector,staggered)
makeBoundaryConditionBase(vertexVector,staggered)
makeBoundaryConditionBase(tensor,staggered)
makeBoundaryConditionBase(sphericalTensor,staggered)
makeBoundaryConditionBase(symmTensor,staggered)
makeBoundaryConditionBase(diagTensor,staggered)
makeBoundaryConditionBase(stencil,staggered)
makeBoundaryConditionBase(diagStencil,staggered)

}

}

}
