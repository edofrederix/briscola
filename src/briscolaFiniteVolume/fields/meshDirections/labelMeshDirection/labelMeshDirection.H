#ifndef labelMeshDirection_H
#define labelMeshDirection_H

#include "label.H"
#include "meshDirection.H"

#define TEMPLATE template<class MeshType>
#include "meshDirectionFunctionsM.H"

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
    #include "labelMeshDirection.C"
#endif

#endif
