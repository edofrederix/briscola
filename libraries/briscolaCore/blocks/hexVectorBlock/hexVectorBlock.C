#include "hexVectorBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

BINARY_TYPE_OPERATOR(hexVector, hexVector, vector, +, add)
BINARY_TYPE_OPERATOR(hexVector, hexVector, hexVector, +, add)
BINARY_TYPE_OPERATOR(hexVector, vector, hexVector, -, subtract)
BINARY_TYPE_OPERATOR(hexVector, hexVector, hexVector, -, subtract)

UNARY_FUNCTION(hexScalar, hexVector, magSqr)
UNARY_FUNCTION(hexScalar, hexVector, mag)
UNARY_FUNCTION(hexScalar, hexVector, cmptMax)
UNARY_FUNCTION(hexScalar, hexVector, cmptMin)
UNARY_FUNCTION(hexScalar, hexVector, cmptSum)
UNARY_FUNCTION(hexScalar, hexVector, cmptAv)
UNARY_FUNCTION(hexScalar, hexVector, cmptProduct)
UNARY_FUNCTION(hexVector, hexVector, cmptSqr)
UNARY_FUNCTION(hexVector, hexVector, cmptMag)

BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptMultiply)
BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptPow)
BINARY_FUNCTION(hexVector, hexVector, hexVector, cmptDivide)
BINARY_FUNCTION(hexVector, hexVector, hexVector, max)
BINARY_FUNCTION(hexVector, hexVector, hexVector, min)
BINARY_FUNCTION(hexVector, hexVector, hexVector, minMod)

}

}

#include "undefBlockFunctionsM.H"