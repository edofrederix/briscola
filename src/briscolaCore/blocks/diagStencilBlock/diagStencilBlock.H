#ifndef diagStencilBlock_H
#define diagStencilBlock_H

#include "block.H"
#include "diagStencil.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "diagStencilBlock.C"
#endif

#endif
