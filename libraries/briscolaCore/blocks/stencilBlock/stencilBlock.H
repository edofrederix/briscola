#ifndef stencilBlock_H
#define stencilBlock_H

#include "scalarBlock.H"
#include "stencil.H"
#include "diagStencil.H"
#include "symmStencil.H"
#include "scalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<stencil> stencilBlock;

BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, stencil, symmStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)
BINARY_TYPE_OPERATOR(stencil, symmStencil, stencil, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#endif
