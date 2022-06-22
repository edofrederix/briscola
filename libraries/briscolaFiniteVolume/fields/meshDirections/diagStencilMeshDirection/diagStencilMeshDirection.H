#ifndef diagStencilMeshDirection_H
#define diagStencilMeshDirection_H

#include "diagStencil.H"
#include "diagStencilBlock.H"
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

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "diagStencilMeshDirection.C"
#endif

#endif
