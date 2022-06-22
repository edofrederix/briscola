#include "hexScalarBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, +, add)
BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, -, subtract)

BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, *, multiply)
BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(hexScalar, hexScalar, hexScalar, /, divide)

}

}
