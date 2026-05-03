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
makeBoundaryConditionBase(vertexScalar,staggered)
makeBoundaryConditionBase(vector,staggered)
makeBoundaryConditionBase(faceVector,staggered)
makeBoundaryConditionBase(vertexVector,staggered)
makeBoundaryConditionBase(tensor,staggered)
makeBoundaryConditionBase(symmTensor,staggered)
makeBoundaryConditionBase(stencil,staggered)
makeBoundaryConditionBase(diagStencil,staggered)

}

}

}
