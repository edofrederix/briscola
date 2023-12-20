#ifndef lowerFaceVectorMeshField_H
#define lowerFaceVectorMeshField_H

#include "lowerFaceVector.H"
#include "scalarMeshField.H"
#include "lowerFaceScalarMeshField.H"
#include "vectorMeshField.H"
#include "tensorMeshField.H"
#include "sphericalTensorMeshField.H"
#include "symmTensorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, lowerFaceVector, +, add)
BINARY_TYPE_OPERATOR(lowerFaceVector, faceVector, lowerFaceVector, +, add)
BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, faceVector, +, add)
BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, vector, +, add)
BINARY_TYPE_OPERATOR(lowerFaceVector, vector, lowerFaceVector, +, add)

BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, lowerFaceVector, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceVector, faceVector, lowerFaceVector, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, faceVector, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceVector, lowerFaceVector, vector, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceVector, vector, lowerFaceVector, -, subtract)

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
BINARY_FUNCTION(lowerFaceVector, scalar, lowerFaceVector, dot)
BINARY_FUNCTION(lowerFaceVector, lowerFaceVector, scalar, dot)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "lowerFaceVectorMeshField.C"
#endif

#endif
