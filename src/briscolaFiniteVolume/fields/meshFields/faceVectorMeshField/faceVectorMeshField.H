#ifndef faceVectorMeshField_H
#define faceVectorMeshField_H

#include "faceVector.H"
#include "scalarMeshField.H"
#include "faceScalarMeshField.H"
#include "vectorMeshField.H"
#include "tensorMeshField.H"
#include "symmTensorMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(faceVector, faceVector, faceVector, +, add)
BINARY_TYPE_OPERATOR(faceVector, faceVector, vector, +, add)
BINARY_TYPE_OPERATOR(faceVector, vector, faceVector, +, add)

BINARY_TYPE_OPERATOR(faceVector, faceVector, faceVector, -, subtract)
BINARY_TYPE_OPERATOR(faceVector, faceVector, vector, -, subtract)
BINARY_TYPE_OPERATOR(faceVector, vector, faceVector, -, subtract)

BINARY_OPERATOR(faceVector, faceVector, faceScalar, /, divide)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, max)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, max)

BINARY_TYPE_FUNCTION(faceVector, faceVector, vector, min)
BINARY_TYPE_FUNCTION(faceVector, vector, faceVector, min)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "faceVectorMeshField.C"
#endif

#endif
