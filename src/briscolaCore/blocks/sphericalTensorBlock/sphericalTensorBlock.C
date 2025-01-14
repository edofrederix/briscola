#include "sphericalTensorBlock.H"

#define TEMPLATE
#include "blockFunctionsM.C"

namespace Foam
{

namespace briscola
{

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)

}

}

#include "undefBlockFunctionsM.H"
