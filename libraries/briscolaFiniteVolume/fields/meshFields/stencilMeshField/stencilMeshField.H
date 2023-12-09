#ifndef stencilMeshField_H
#define stencilMeshField_H

#include "scalarMeshField.H"
#include "vectorMeshField.H"
#include "diagStencilMeshField.H"

#define TEMPLATE template<class MeshType>
#include "meshFieldFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(stencil, stencil, symmStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)

BINARY_TYPE_OPERATOR(stencil, symmStencil, stencil, -, subtract)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "stencilMeshField.C"
#endif

#endif
