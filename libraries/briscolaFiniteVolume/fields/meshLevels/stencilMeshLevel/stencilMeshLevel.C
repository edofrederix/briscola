#include "stencilMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(stencil, stencil, symmStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)

BINARY_TYPE_OPERATOR(stencil, symmStencil, stencil, -, subtract)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"
