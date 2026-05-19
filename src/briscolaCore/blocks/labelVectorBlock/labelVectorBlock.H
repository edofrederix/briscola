#ifndef labelVectorBlock_H
#define labelVectorBlock_H

#include "block.H"
#include "labelVector.H"

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
    #include "labelVectorBlock.C"
#endif

#endif
