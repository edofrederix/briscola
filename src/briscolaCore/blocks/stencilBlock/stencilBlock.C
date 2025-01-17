#include "stencilBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

#include "undefBlockFunctionsM.H"
