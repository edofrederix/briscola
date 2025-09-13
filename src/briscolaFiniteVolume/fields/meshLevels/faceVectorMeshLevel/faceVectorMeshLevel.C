#include "faceVectorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_OPERATOR(faceVector, faceVector, faceScalar, /, divide)

UNARY_FUNCTION(faceScalar, faceVector, magSqr)
UNARY_FUNCTION(faceScalar, faceVector, mag)
UNARY_FUNCTION(faceScalar, faceVector, cmptMax)
UNARY_FUNCTION(faceScalar, faceVector, cmptMin)
UNARY_FUNCTION(faceScalar, faceVector, cmptSum)
UNARY_FUNCTION(faceScalar, faceVector, cmptAv)
UNARY_FUNCTION(faceScalar, faceVector, cmptProduct)
UNARY_FUNCTION(faceVector, faceVector, cmptSqr)
UNARY_FUNCTION(faceVector, faceVector, cmptMag)

BINARY_FUNCTION(faceVector, faceVector, faceVector, cmptMultiply)
BINARY_FUNCTION(faceVector, faceVector, faceVector, cmptPow)
BINARY_FUNCTION(faceVector, faceVector, faceVector, cmptDivide)
BINARY_FUNCTION(faceVector, faceVector, faceVector, max)
BINARY_FUNCTION(faceVector, faceVector, faceVector, min)
BINARY_FUNCTION(faceVector, faceVector, faceVector, minMod)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, max)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, max)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, min)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, min)

}

}

}

#include "undefBlockFunctionsM.H"

