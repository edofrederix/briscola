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

UNARY_FUNCTION(vertexScalar, vertexScalar, sqr)
UNARY_FUNCTION(vertexScalar, vertexScalar, sqrt)
UNARY_FUNCTION(vertexScalar, vertexScalar, cbrt)
UNARY_FUNCTION(vertexScalar, vertexScalar, sign)
UNARY_FUNCTION(vertexScalar, vertexScalar, pos)
UNARY_FUNCTION(vertexScalar, vertexScalar, pos0)
UNARY_FUNCTION(vertexScalar, vertexScalar, neg)
UNARY_FUNCTION(vertexScalar, vertexScalar, neg0)
UNARY_FUNCTION(vertexScalar, vertexScalar, exp)
UNARY_FUNCTION(vertexScalar, vertexScalar, log)
UNARY_FUNCTION(vertexScalar, vertexScalar, sin)
UNARY_FUNCTION(vertexScalar, vertexScalar, cos)
UNARY_FUNCTION(vertexScalar, vertexScalar, tan)
UNARY_FUNCTION(vertexScalar, vertexScalar, asin)
UNARY_FUNCTION(vertexScalar, vertexScalar, acos)
UNARY_FUNCTION(vertexScalar, vertexScalar, atan)
UNARY_FUNCTION(vertexScalar, vertexScalar, erf)

}

}

}

#include "undefBlockFunctionsM.H"
