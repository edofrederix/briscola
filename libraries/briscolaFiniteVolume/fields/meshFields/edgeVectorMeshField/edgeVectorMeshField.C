#include "edgeVectorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(edgeVector, edgeVector, edgeVector, +, add)
BINARY_TYPE_OPERATOR(edgeVector, vector, edgeVector, +, add)

BINARY_TYPE_OPERATOR(edgeVector, edgeVector, edgeVector, -, subtract)
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
BINARY_FUNCTION(edgeVector, scalar, edgeVector, dot)
BINARY_FUNCTION(edgeVector, edgeVector, scalar, dot)

}

}

}

#include "undefBlockFunctionsM.H"
