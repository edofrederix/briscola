#ifndef tensorBlock_H
#define tensorBlock_H

#include "scalarBlock.H"
#include "vectorBlock.H"
#include "symmTensorBlock.H"
#include "tensor.H"

#define TEMPLATE
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

typedef block<tensor> tensorBlock;

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

#endif
