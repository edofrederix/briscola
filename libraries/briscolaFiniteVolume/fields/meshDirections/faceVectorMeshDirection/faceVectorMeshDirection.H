#ifndef faceVectorMeshDirection_H
#define faceVectorMeshDirection_H

#include "faceVector.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"
#include "tensorMeshDirection.H"
#include "sphericalTensorMeshDirection.H"
#include "symmTensorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

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

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "faceVectorMeshDirection.C"
#endif

#endif
