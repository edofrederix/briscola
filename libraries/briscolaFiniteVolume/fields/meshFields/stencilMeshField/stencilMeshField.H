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

}

}

}

#include "undefBlockFunctionsM.H"

#ifdef NoRepository
    #include "stencilMeshField.C"
#endif

#endif
