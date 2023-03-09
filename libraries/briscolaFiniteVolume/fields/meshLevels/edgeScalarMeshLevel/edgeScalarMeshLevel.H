#ifndef edgeScalarMeshLevel_H
#define edgeScalarMeshLevel_H

#include "edgeScalar.H"
#include "meshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "edgeScalarMeshLevel.C"
#endif

#endif
