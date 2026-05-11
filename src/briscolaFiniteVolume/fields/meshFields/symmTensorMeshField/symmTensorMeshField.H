#ifndef symmTensorMeshField_H
#define symmTensorMeshField_H

#include "meshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

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

#ifdef NoRepository
    #include "symmTensorMeshField.C"
#endif

#endif
