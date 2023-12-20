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
BINARY_TYPE_OPERATOR(faceScalar, faceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, +, add)

BINARY_TYPE_OPERATOR(faceScalar, faceScalar, faceScalar, -, subtract)
BINARY_TYPE_OPERATOR(faceScalar, faceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, -, subtract)

BINARY_OPERATOR(faceScalar, faceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(faceScalar, scalar, faceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceScalarMeshDirection.C"
#endif

#endif
