#include "lowerFaceScalarBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, *, multiply)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

}

}

#include "undefBlockFunctionsM.H"