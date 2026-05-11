#include "symmTensorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.C"

namespace Foam
{

namespace briscola
{

namespace fv
{

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(scalar, symmTensor, det)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)

}

}

}

#include "undefBlockFunctionsM.H"
