#ifndef vertexVectorMeshDirection_H
#define vertexVectorMeshDirection_H

#include "vertexVector.H"
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

BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"


#ifdef NoRepository
    #include "vertexVectorMeshDirection.C"
#endif

#endif
