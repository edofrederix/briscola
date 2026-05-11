#ifndef faceVectorMeshLevel_H
#define faceVectorMeshLevel_H

#include "faceVector.H"
#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "vectorMeshLevel.H"
#include "tensorMeshLevel.H"
#include "symmTensorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "faceVectorMeshLevel.C"
#endif

#endif
