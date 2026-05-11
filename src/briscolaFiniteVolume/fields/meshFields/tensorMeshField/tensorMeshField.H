#ifndef tensorMeshField_H
#define tensorMeshField_H

#include "symmTensorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

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
    #include "tensorMeshField.C"
#endif

#endif
