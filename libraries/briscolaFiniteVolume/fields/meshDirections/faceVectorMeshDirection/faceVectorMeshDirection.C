#include "faceVectorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(faceVector, vector, faceVector, +, add)
BINARY_TYPE_OPERATOR(faceVector, vector, faceVector, -, subtract)

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

}

}

}

#include "undefBlockFunctionsM.H"

