#ifndef edgeVectorMeshLevel_H
#define edgeVectorMeshLevel_H

#include "edgeVector.H"
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

BINARY_TYPE_OPERATOR(edgeVector, vector, edgeVector, +, add)
BINARY_TYPE_OPERATOR(edgeVector, vector, edgeVector, -, subtract)

UNARY_FUNCTION(edgeScalar, edgeVector, magSqr)
UNARY_FUNCTION(edgeScalar, edgeVector, mag)
UNARY_FUNCTION(edgeScalar, edgeVector, cmptMax)
UNARY_FUNCTION(edgeScalar, edgeVector, cmptMin)
UNARY_FUNCTION(edgeScalar, edgeVector, cmptSum)
UNARY_FUNCTION(edgeScalar, edgeVector, cmptAv)
UNARY_FUNCTION(edgeScalar, edgeVector, cmptProduct)
UNARY_FUNCTION(edgeVector, edgeVector, cmptSqr)
UNARY_FUNCTION(edgeVector, edgeVector, cmptMag)

BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, cmptMultiply)
BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, cmptPow)
BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, cmptDivide)
BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, max)
BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, min)
BINARY_FUNCTION(edgeVector, edgeVector, edgeVector, minMod)

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "edgeVectorMeshLevel.C"
#endif

#endif
