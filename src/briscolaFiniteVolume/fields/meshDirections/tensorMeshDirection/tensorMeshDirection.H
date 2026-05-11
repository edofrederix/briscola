#ifndef tensorMeshDirection_H
#define tensorMeshDirection_H

#include "tensor.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"
#include "symmTensorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
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

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "tensorMeshDirection.C"
#endif

#endif
