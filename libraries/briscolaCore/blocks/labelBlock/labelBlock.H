#ifndef labelBlock_H
#define labelBlock_H

#include "label.H"
#include "block.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<label> labelBlock;

UNARY_FUNCTION(label, label, mag)
BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)

}

}

#include "undefBlockFunctionsM.H"

#endif
