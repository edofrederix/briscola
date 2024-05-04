#include "lowerFaceScalarMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"
