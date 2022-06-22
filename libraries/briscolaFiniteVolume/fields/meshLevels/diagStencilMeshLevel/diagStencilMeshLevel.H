#ifndef diagStencilMeshLevel_H
#define diagStencilMeshLevel_H

#include "diagStencil.H"
#include "diagStencilBlock.H"
#include "vector.H"

#include "meshLevel.H"
#include "scalarMeshLevel.H"
#include "vectorMeshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

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
    #include "diagStencilMeshLevel.C"
#endif

#endif
