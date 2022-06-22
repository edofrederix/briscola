#ifndef stencilMeshDirection_H
#define stencilMeshDirection_H

#include "stencil.H"
#include "vector.H"
#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"
#include "symmStencilMeshDirection.H"
#include "diagStencilMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(stencil, stencil, symmStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, symmStencil, stencil, -, subtract)

BINARY_TYPE_OPERATOR(stencil, stencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(stencil, diagStencil, stencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "stencilMeshDirection.C"
#endif

#endif
