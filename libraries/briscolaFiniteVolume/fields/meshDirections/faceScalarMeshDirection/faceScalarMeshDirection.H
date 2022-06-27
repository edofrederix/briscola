#ifndef faceScalarMeshDirection_H
#define faceScalarMeshDirection_H

#include "faceScalar.H"
#include "meshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(faceScalar, faceScalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, faceScalar, faceScalar, -, subtract)

BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, *, multiply)
BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(faceScalar, faceScalar, faceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceScalarMeshDirection.C"
#endif

#endif
