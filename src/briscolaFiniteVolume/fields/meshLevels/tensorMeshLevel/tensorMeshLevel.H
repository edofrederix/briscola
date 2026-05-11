#ifndef tensorMeshLevel_H
#define tensorMeshLevel_H

#include "tensor.H"
#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "vectorMeshLevel.H"
#include "symmTensorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "tensorMeshLevel.C"
#endif

#endif
