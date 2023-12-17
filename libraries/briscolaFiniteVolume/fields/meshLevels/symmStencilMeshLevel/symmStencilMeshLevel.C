#include "symmStencilMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(symmStencil, symmStencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, stencil, symmStencil, +, add)

BINARY_TYPE_OPERATOR(symmStencil, diagStencil, symmStencil, -, subtract)
BINARY_TYPE_OPERATOR(symmStencil, symmStencil, stencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"
