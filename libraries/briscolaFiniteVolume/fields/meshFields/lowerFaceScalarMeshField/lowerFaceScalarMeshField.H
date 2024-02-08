#ifndef lowerFaceScalarMeshField_H
#define lowerFaceScalarMeshField_H

#include "meshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "lowerFaceScalarMeshField.C"
#endif

#endif
