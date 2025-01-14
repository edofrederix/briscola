#ifndef sphericalTensorMeshDirection_H
#define sphericalTensorMeshDirection_H

#include "sphericalTensor.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "sphericalTensorMeshDirection.C"
#endif

#endif
