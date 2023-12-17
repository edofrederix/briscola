#ifndef symmStencilBlock_H
#define symmStencilBlock_H

#include "scalarBlock.H"
#include "symmStencil.H"
#include "diagStencil.H"
#include "scalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<symmStencil> symmStencilBlock;

BINARY_TYPE_OPERATOR(symmStencil, symmStencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, stencil, symmStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, diagStencil, symmStencil, -, subtract)
BINARY_TYPE_OPERATOR(symmStencil, symmStencil, stencil, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#endif
