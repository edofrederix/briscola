#ifndef vertexVectorMeshLevel_H
#define vertexVectorMeshLevel_H

#include "vertexVector.H"
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

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "vertexVectorMeshLevel.C"
#endif

#endif
