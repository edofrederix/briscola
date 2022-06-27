#include "faceScalarBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, -, subtract)

BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, *, multiply)
BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(faceScalar, faceScalar, faceScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"