#ifndef hexScalarMeshField_H
#define hexScalarMeshField_H

#include "meshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, +, add)
BINARY_TYPE_OPERATOR(hexScalar, hexScalar, hexScalar, -, subtract)

BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, *, multiply)
BINARY_OPERATOR(hexScalar, hexScalar, hexScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(hexScalar, hexScalar, hexScalar, /, divide)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "hexScalarMeshField.C"
#endif

#endif
