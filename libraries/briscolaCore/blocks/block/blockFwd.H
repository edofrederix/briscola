#ifndef blockFwd_H
#define blockFwd_H

#include "label.H"
#include "scalar.H"
#include "hexScalar.H"
#include "vector.H"
#include "hexVector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "diagTensor.H"
#include "symmStencil.H"
#include "diagStencil.H"
#include "stencil.H"

namespace Foam
{

namespace briscola
{

template<class Type> class block;

typedef block<label> labelBlock;
typedef block<scalar> scalarBlock;
typedef block<hexScalar> hexScalarBlock;
typedef block<vector> vectorBlock;
typedef block<hexVector> hexVectorBlock;
typedef block<tensor> tensorBlock;
typedef block<symmTensor> symmTensorBlock;
typedef block<sphericalTensor> sphericalTensorBlock;
typedef block<diagTensor> diagTensorBlock;
typedef block<symmStencil> symmStencilBlock;
typedef block<diagStencil> diagStencilBlock;
typedef block<stencil> stencilBlock;

}

}

#endif
