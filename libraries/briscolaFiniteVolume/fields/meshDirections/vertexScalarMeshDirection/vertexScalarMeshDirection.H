#ifndef vertexScalarMeshDirection_H
#define vertexScalarMeshDirection_H

#include "vertexScalar.H"
#include "meshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, +, add)
BINARY_TYPE_OPERATOR(vertexScalar, vertexScalar, vertexScalar, -, subtract)

BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, *, multiply)
BINARY_OPERATOR(vertexScalar, vertexScalar, vertexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(vertexScalar, vertexScalar, vertexScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "vertexScalarMeshDirection.C"
#endif

#endif
