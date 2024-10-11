#ifndef edgeScalarBlock_H
#define edgeScalarBlock_H

#include "block.H"
#include "edgeScalar.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<edgeScalar> edgeScalarBlock;

BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, +, add)
BINARY_TYPE_OPERATOR(edgeScalar, scalar, edgeScalar, -, subtract)

BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, *, multiply)
BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(edgeScalar, edgeScalar, edgeScalar, /, divide)

UNARY_FUNCTION(edgeScalar, edgeScalar, pow3)
UNARY_FUNCTION(edgeScalar, edgeScalar, pow4)
UNARY_FUNCTION(edgeScalar, edgeScalar, pow5)
UNARY_FUNCTION(edgeScalar, edgeScalar, pow6)
UNARY_FUNCTION(edgeScalar, edgeScalar, pow025)
UNARY_FUNCTION(edgeScalar, edgeScalar, sqrt)
UNARY_FUNCTION(edgeScalar, edgeScalar, cbrt)
UNARY_FUNCTION(edgeScalar, edgeScalar, sign)
UNARY_FUNCTION(edgeScalar, edgeScalar, pos)
UNARY_FUNCTION(edgeScalar, edgeScalar, pos0)
UNARY_FUNCTION(edgeScalar, edgeScalar, neg)
UNARY_FUNCTION(edgeScalar, edgeScalar, neg0)
UNARY_FUNCTION(edgeScalar, edgeScalar, posPart)
UNARY_FUNCTION(edgeScalar, edgeScalar, negPart)
UNARY_FUNCTION(edgeScalar, edgeScalar, exp)
UNARY_FUNCTION(edgeScalar, edgeScalar, log)
UNARY_FUNCTION(edgeScalar, edgeScalar, log10)
UNARY_FUNCTION(edgeScalar, edgeScalar, sin)
UNARY_FUNCTION(edgeScalar, edgeScalar, cos)
UNARY_FUNCTION(edgeScalar, edgeScalar, tan)
UNARY_FUNCTION(edgeScalar, edgeScalar, asin)
UNARY_FUNCTION(edgeScalar, edgeScalar, acos)
UNARY_FUNCTION(edgeScalar, edgeScalar, atan)
UNARY_FUNCTION(edgeScalar, edgeScalar, sinh)
UNARY_FUNCTION(edgeScalar, edgeScalar, cosh)
UNARY_FUNCTION(edgeScalar, edgeScalar, tanh)
UNARY_FUNCTION(edgeScalar, edgeScalar, asinh)
UNARY_FUNCTION(edgeScalar, edgeScalar, acosh)
UNARY_FUNCTION(edgeScalar, edgeScalar, atanh)
UNARY_FUNCTION(edgeScalar, edgeScalar, erf)
UNARY_FUNCTION(edgeScalar, edgeScalar, erfc)
UNARY_FUNCTION(edgeScalar, edgeScalar, lgamma)
UNARY_FUNCTION(edgeScalar, edgeScalar, j0)
UNARY_FUNCTION(edgeScalar, edgeScalar, j1)
UNARY_FUNCTION(edgeScalar, edgeScalar, y0)
UNARY_FUNCTION(edgeScalar, edgeScalar, y1)

}

}

#include "undefBlockFunctionsM.H"

#endif
