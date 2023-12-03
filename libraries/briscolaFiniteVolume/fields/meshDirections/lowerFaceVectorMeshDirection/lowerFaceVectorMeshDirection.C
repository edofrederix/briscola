#include "lowerFaceVectorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceVector, vector, lowerFaceVector, +, add)
BINARY_TYPE_OPERATOR(lowerFaceVector, vector, lowerFaceVector, -, subtract)

UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, magSqr)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, mag)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptMax)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptMin)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptSum)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptAv)
UNARY_FUNCTION(lowerFaceScalar, lowerFaceVector, cmptProduct)
UNARY_FUNCTION(lowerFaceVector, lowerFaceVector, cmptSqr)
UNARY_FUNCTION(lowerFaceVector, lowerFaceVector, cmptMag)

BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptMultiply)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptPow)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, cmptDivide)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, max)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, min)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, lowerFaceVector, minMod)

}

}

}

#include "undefBlockFunctionsM.H"

