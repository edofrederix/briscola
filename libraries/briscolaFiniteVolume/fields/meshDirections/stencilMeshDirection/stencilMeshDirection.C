#include "stencilMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"
