#ifndef edgeScalarMeshField_H
#define edgeScalarMeshField_H

#include "meshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(edgeScalar, edgeScalar, edgeScalar, +, add)
BINARY_TYPE_OPERATOR(edgeScalar, edgeScalar, edgeScalar, -, subtract)

BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, *, multiply)
BINARY_OPERATOR(edgeScalar, edgeScalar, edgeScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(edgeScalar, edgeScalar, edgeScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "edgeScalarMeshField.C"
#endif

#endif
