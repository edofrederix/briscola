#ifndef stencilMeshLevel_H
#define stencilMeshLevel_H

#include "stencil.H"
#include "vector.H"
#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "vectorMeshLevel.H"
#include "diagStencilMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "stencilMeshLevel.C"
#endif

#endif
