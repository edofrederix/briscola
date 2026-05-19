#ifndef stencilBlock_H
#define stencilBlock_H

#include "scalarBlock.H"
#include "stencil.H"
#include "diagStencil.H"
#include "scalar.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "stencilBlock.C"
#endif

#endif
