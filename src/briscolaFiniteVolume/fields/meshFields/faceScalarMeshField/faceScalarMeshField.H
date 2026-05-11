#ifndef faceScalarMeshField_H
#define faceScalarMeshField_H

#include "meshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

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

BINARY_TYPE_FUNCTION_FS(faceScalar, faceScalar, scalar, pow)

UNARY_FUNCTION(faceScalar, faceScalar, sqr)
UNARY_FUNCTION(faceScalar, faceScalar, sqrt)
UNARY_FUNCTION(faceScalar, faceScalar, cbrt)
UNARY_FUNCTION(faceScalar, faceScalar, sign)
UNARY_FUNCTION(faceScalar, faceScalar, pos)
UNARY_FUNCTION(faceScalar, faceScalar, pos0)
UNARY_FUNCTION(faceScalar, faceScalar, neg)
UNARY_FUNCTION(faceScalar, faceScalar, neg0)
UNARY_FUNCTION(faceScalar, faceScalar, exp)
UNARY_FUNCTION(faceScalar, faceScalar, log)
UNARY_FUNCTION(faceScalar, faceScalar, sin)
UNARY_FUNCTION(faceScalar, faceScalar, cos)
UNARY_FUNCTION(faceScalar, faceScalar, tan)
UNARY_FUNCTION(faceScalar, faceScalar, asin)
UNARY_FUNCTION(faceScalar, faceScalar, acos)
UNARY_FUNCTION(faceScalar, faceScalar, atan)
UNARY_FUNCTION(faceScalar, faceScalar, erf)

BINARY_TYPE_FUNCTION(faceScalar, faceScalar, scalar, max)
BINARY_TYPE_FUNCTION(faceScalar, scalar, faceScalar, max)

BINARY_TYPE_FUNCTION(faceScalar, faceScalar, scalar, min)
BINARY_TYPE_FUNCTION(faceScalar, scalar, faceScalar, min)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceScalarMeshField.C"
#endif

#endif
