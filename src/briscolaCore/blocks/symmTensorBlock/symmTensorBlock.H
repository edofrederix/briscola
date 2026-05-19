#ifndef symmTensorBlock_H
#define symmTensorBlock_H

#include "block.H"
#include "symmTensor.H"
#include "tensor.H"

#define TEMPLATE template<int P>
#include "blockFunctionsM.H"

namespace Foam
{

namespace briscola
{

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(scalar, symmTensor, det)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "symmTensorBlock.C"
#endif

#endif
