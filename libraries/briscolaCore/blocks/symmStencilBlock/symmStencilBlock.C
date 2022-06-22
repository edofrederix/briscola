#include "symmStencilBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(symmStencil, symmStencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, diagStencil, symmStencil, -, subtract)

}

}

#include "undefBlockFunctionsM.H"