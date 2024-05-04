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

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, +, add)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, +, add)

BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, lowerFaceScalar, scalar, -, subtract)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, -, subtract)

BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, faceScalar, lowerFaceScalar, /, divide)
BINARY_OPERATOR(lowerFaceScalar, lowerFaceScalar, faceScalar, /, divide)
BINARY_TYPE_OPERATOR(lowerFaceScalar, scalar, lowerFaceScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "lowerFaceScalarMeshField.C"
#endif

#endif
