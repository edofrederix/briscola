#ifndef vertexVectorMeshField_H
#define vertexVectorMeshField_H

#include "vertexVector.H"
#include "scalarMeshField.H"
#include "vertexScalarMeshField.H"
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

BINARY_TYPE_OPERATOR(vertexVector, vertexVector, vertexVector, +, add)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, +, add)

BINARY_TYPE_OPERATOR(vertexVector, vertexVector, vertexVector, -, subtract)
BINARY_TYPE_OPERATOR(vertexVector, vector, vertexVector, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "vertexVectorMeshField.C"
#endif

#endif
