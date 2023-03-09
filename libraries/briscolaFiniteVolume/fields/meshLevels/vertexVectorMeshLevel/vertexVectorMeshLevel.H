#ifndef vertexVectorMeshLevel_H
#define vertexVectorMeshLevel_H

#include "vertexVector.H"
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

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

UNARY_FUNCTION(vertexScalar, vertexVector, magSqr)
UNARY_FUNCTION(vertexScalar, vertexVector, mag)
UNARY_FUNCTION(vertexScalar, vertexVector, cmptMax)
UNARY_FUNCTION(vertexScalar, vertexVector, cmptMin)
UNARY_FUNCTION(vertexScalar, vertexVector, cmptSum)
UNARY_FUNCTION(vertexScalar, vertexVector, cmptAv)
UNARY_FUNCTION(vertexScalar, vertexVector, cmptProduct)
UNARY_FUNCTION(vertexVector, vertexVector, cmptSqr)
UNARY_FUNCTION(vertexVector, vertexVector, cmptMag)

BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, cmptMultiply)
BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, cmptPow)
BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, cmptDivide)
BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, max)
BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, min)
BINARY_FUNCTION(vertexVector, vertexVector, vertexVector, minMod)

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "vertexVectorMeshLevel.C"
#endif

#endif
