#ifndef vectorMeshLevel_H
#define vectorMeshLevel_H

#include "vector.H"
#include "meshLevel.H"

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
    #include "vectorMeshLevel.C"
#endif

#endif
