#include "vertexVectorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(vertexVector, vertexVector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)

BINARY_TYPE_OPERATOR(vertexVector, vertexVector, vertexVector, -, subtract)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"
