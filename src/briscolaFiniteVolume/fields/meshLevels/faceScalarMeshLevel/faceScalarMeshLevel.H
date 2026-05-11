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

BINARY_TYPE_FUNCTION_FS(faceScalar, faceScalar, scalar, pow)

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

BINARY_TYPE_FUNCTION(faceScalar, faceScalar, scalar, max)
BINARY_TYPE_FUNCTION(faceScalar, scalar, faceScalar, max)

BINARY_TYPE_FUNCTION(faceScalar, faceScalar, scalar, min)
BINARY_TYPE_FUNCTION(faceScalar, scalar, faceScalar, min)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceScalarMeshLevel.C"
#endif

#endif
