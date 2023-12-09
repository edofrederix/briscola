#ifndef lowerFaceVectorMeshLevel_H
#define lowerFaceVectorMeshLevel_H

#include "lowerFaceVector.H"
#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "vectorMeshLevel.H"
#include "tensorMeshLevel.H"
#include "sphericalTensorMeshLevel.H"
#include "symmTensorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_OPERATOR(lowerFaceVector, lowerFaceVector, lowerFaceScalar, /, divide)

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


#ifdef NoRepository
    #include "lowerFaceVectorMeshLevel.C"
#endif

#endif
