#include "edgeScalarBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, +, add)
BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, -, subtract)

BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, *, multiply)
BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(edgeScalar, edgeScalar, edgeScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"