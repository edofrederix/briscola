#ifndef faceScalarMeshLevel_H
#define faceScalarMeshLevel_H

#include "faceScalar.H"
#include "meshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(faceScalar, faceScalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, faceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, +, add)

BINARY_TYPE_OPERATOR(faceScalar, faceScalar, faceScalar, -, subtract)
BINARY_TYPE_OPERATOR(faceScalar, faceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, -, subtract)

BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, /, divide)

UNARY_FUNCTION(faceScalar, faceScalar, pow3)
UNARY_FUNCTION(faceScalar, faceScalar, pow4)
UNARY_FUNCTION(faceScalar, faceScalar, pow5)
UNARY_FUNCTION(faceScalar, faceScalar, pow6)
UNARY_FUNCTION(faceScalar, faceScalar, pow025)
UNARY_FUNCTION(faceScalar, faceScalar, sqrt)
UNARY_FUNCTION(faceScalar, faceScalar, cbrt)
UNARY_FUNCTION(faceScalar, faceScalar, sign)
UNARY_FUNCTION(faceScalar, faceScalar, pos)
UNARY_FUNCTION(faceScalar, faceScalar, pos0)
UNARY_FUNCTION(faceScalar, faceScalar, neg)
UNARY_FUNCTION(faceScalar, faceScalar, neg0)
UNARY_FUNCTION(faceScalar, faceScalar, posPart)
UNARY_FUNCTION(faceScalar, faceScalar, negPart)
UNARY_FUNCTION(faceScalar, faceScalar, exp)
UNARY_FUNCTION(faceScalar, faceScalar, log)
UNARY_FUNCTION(faceScalar, faceScalar, log10)
UNARY_FUNCTION(faceScalar, faceScalar, sin)
UNARY_FUNCTION(faceScalar, faceScalar, cos)
UNARY_FUNCTION(faceScalar, faceScalar, tan)
UNARY_FUNCTION(faceScalar, faceScalar, asin)
UNARY_FUNCTION(faceScalar, faceScalar, acos)
UNARY_FUNCTION(faceScalar, faceScalar, atan)
UNARY_FUNCTION(faceScalar, faceScalar, sinh)
UNARY_FUNCTION(faceScalar, faceScalar, cosh)
UNARY_FUNCTION(faceScalar, faceScalar, tanh)
UNARY_FUNCTION(faceScalar, faceScalar, asinh)
UNARY_FUNCTION(faceScalar, faceScalar, acosh)
UNARY_FUNCTION(faceScalar, faceScalar, atanh)
UNARY_FUNCTION(faceScalar, faceScalar, erf)
UNARY_FUNCTION(faceScalar, faceScalar, erfc)
UNARY_FUNCTION(faceScalar, faceScalar, lgamma)
UNARY_FUNCTION(faceScalar, faceScalar, j0)
UNARY_FUNCTION(faceScalar, faceScalar, j1)
UNARY_FUNCTION(faceScalar, faceScalar, y0)
UNARY_FUNCTION(faceScalar, faceScalar, y1)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceScalarMeshLevel.C"
#endif

#endif
