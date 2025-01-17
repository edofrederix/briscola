#include "labelMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

UNARY_FUNCTION(label, label, mag)
BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"
