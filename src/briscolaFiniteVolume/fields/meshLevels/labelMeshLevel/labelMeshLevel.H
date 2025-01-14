#ifndef labelMeshLevel_H
#define labelMeshLevel_H

#include "label.H"
#include "meshLevel.H"

#define TEMPLATE template<class MeshType>
#include "meshLevelFunctionsM.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

UNARY_FUNCTION(label, label, mag)
BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "labelMeshLevel.C"
#endif

#endif
