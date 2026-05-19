#include "tensorBlock.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(tensor, tensor, T)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(scalar, tensor, det)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)

}

}

#include "undefBlockFunctionsM.H"
