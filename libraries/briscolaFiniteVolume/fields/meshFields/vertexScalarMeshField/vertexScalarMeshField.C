#include "vertexScalarMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, +, add)
BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, -, subtract)

BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, *, multiply)
BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(vertexScalar, vertexScalar, vertexScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"
