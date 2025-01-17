#ifndef sphericalTensorBlock_H
#define sphericalTensorBlock_H

#include "scalarBlock.H"
#include "sphericalTensor.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<sphericalTensor> sphericalTensorBlock;

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)

}

}

#include "undefBlockFunctionsM.H"

#endif
