#ifndef symmTensorMeshLevel_H
#define symmTensorMeshLevel_H

#include "symmTensor.H"
#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "tensorMeshLevel.H"
#include "vectorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "symmTensorMeshLevel.C"
#endif

#endif
