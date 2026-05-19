#ifndef labelBlock_H
#define labelBlock_H

#include "label.H"
#include "block.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

UNARY_FUNCTION(label, label, mag)
BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "labelBlock.C"
#endif

#endif
