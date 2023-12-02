#ifndef lowerFaceScalarMeshDirection_H
#define lowerFaceScalarMeshDirection_H

#include "lowerFaceScalar.H"
#include "meshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, *, multiply)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "lowerFaceScalarMeshDirection.C"
#endif

#endif
