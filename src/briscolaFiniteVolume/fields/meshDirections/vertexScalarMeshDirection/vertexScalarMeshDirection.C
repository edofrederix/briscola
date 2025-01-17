#include "vertexScalarMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, +, add)
BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, -, subtract)

BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, *, multiply)
BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(vertexScalar, vertexScalar, vertexScalar, /, divide)

UNARY_FUNCTION(vertexScalar, vertexScalar, pow3)
UNARY_FUNCTION(vertexScalar, vertexScalar, pow4)
UNARY_FUNCTION(vertexScalar, vertexScalar, pow5)
UNARY_FUNCTION(vertexScalar, vertexScalar, pow6)
UNARY_FUNCTION(vertexScalar, vertexScalar, pow025)
UNARY_FUNCTION(vertexScalar, vertexScalar, sqrt)
UNARY_FUNCTION(vertexScalar, vertexScalar, cbrt)
UNARY_FUNCTION(vertexScalar, vertexScalar, sign)
UNARY_FUNCTION(vertexScalar, vertexScalar, pos)
UNARY_FUNCTION(vertexScalar, vertexScalar, pos0)
UNARY_FUNCTION(vertexScalar, vertexScalar, neg)
UNARY_FUNCTION(vertexScalar, vertexScalar, neg0)
UNARY_FUNCTION(vertexScalar, vertexScalar, posPart)
UNARY_FUNCTION(vertexScalar, vertexScalar, negPart)
UNARY_FUNCTION(vertexScalar, vertexScalar, exp)
UNARY_FUNCTION(vertexScalar, vertexScalar, log)
UNARY_FUNCTION(vertexScalar, vertexScalar, log10)
UNARY_FUNCTION(vertexScalar, vertexScalar, sin)
UNARY_FUNCTION(vertexScalar, vertexScalar, cos)
UNARY_FUNCTION(vertexScalar, vertexScalar, tan)
UNARY_FUNCTION(vertexScalar, vertexScalar, asin)
UNARY_FUNCTION(vertexScalar, vertexScalar, acos)
UNARY_FUNCTION(vertexScalar, vertexScalar, atan)
UNARY_FUNCTION(vertexScalar, vertexScalar, sinh)
UNARY_FUNCTION(vertexScalar, vertexScalar, cosh)
UNARY_FUNCTION(vertexScalar, vertexScalar, tanh)
UNARY_FUNCTION(vertexScalar, vertexScalar, asinh)
UNARY_FUNCTION(vertexScalar, vertexScalar, acosh)
UNARY_FUNCTION(vertexScalar, vertexScalar, atanh)
UNARY_FUNCTION(vertexScalar, vertexScalar, erf)
UNARY_FUNCTION(vertexScalar, vertexScalar, erfc)
UNARY_FUNCTION(vertexScalar, vertexScalar, lgamma)
UNARY_FUNCTION(vertexScalar, vertexScalar, j0)
UNARY_FUNCTION(vertexScalar, vertexScalar, j1)
UNARY_FUNCTION(vertexScalar, vertexScalar, y0)
UNARY_FUNCTION(vertexScalar, vertexScalar, y1)

}

}

}

#include "undefBlockFunctionsM.H"
