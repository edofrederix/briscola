#ifndef diagStencilBlock_H
#define diagStencilBlock_H

#include "scalarBlock.H"
#include "diagStencil.H"
#include "scalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<diagStencil> diagStencilBlock;

}

}

#include "undefBlockFunctionsM.H"

#endif
