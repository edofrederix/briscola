#include "lowerFaceScalarMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, *, multiply)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"
