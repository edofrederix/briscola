#ifndef faceVectorMeshDirection_H
#define faceVectorMeshDirection_H

#include "faceVector.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"
#include "tensorMeshDirection.H"
#include "symmTensorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

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
    #include "faceVectorMeshDirection.C"
#endif

#endif
