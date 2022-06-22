#ifndef symmStencilBlock_H
#define symmStencilBlock_H

#include "scalarBlock.H"
#include "symmStencil.H"
#include "scalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<symmStencil> symmStencilBlock;

BINARY_TYPE_OPERATOR(symmStencil, symmStencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, diagStencil, symmStencil, -, subtract)
}

}

#include "undefBlockFunctionsM.H"

#endif
