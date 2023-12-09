#ifndef symmStencilMeshDirection_H
#define symmStencilMeshDirection_H

#include "symmStencil.H"
#include "symmStencilBlock.H"
#include "vector.H"

#include "meshDirection.H"
#include "scalarMeshDirection.H"
#include "vectorMeshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

BINARY_TYPE_OPERATOR(symmStencil, symmStencil, diagStencil, +, add)
BINARY_TYPE_OPERATOR(symmStencil, diagStencil, symmStencil, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "symmStencilMeshDirection.C"
#endif

#endif
