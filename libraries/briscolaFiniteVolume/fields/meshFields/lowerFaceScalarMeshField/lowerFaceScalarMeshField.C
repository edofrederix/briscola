#include "lowerFaceScalarMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, /, divide)

UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow3)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow4)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow5)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow6)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pow025)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sqrt)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cbrt)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sign)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pos)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, pos0)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, neg)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, neg0)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, posPart)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, negPart)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, exp)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, log)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, log10)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sin)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cos)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, tan)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, asin)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, acos)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, atan)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, sinh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, cosh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, tanh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, asinh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, acosh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, atanh)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, erf)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, erfc)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, lgamma)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, j0)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, j1)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, y0)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceScalar, y1)

}

}

}

#include "undefBlockFunctionsM.H"
