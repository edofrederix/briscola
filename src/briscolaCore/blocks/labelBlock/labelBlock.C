#include "labelBlock.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.C"

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
