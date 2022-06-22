#ifndef symmStencilMeshField_H
#define symmStencilMeshField_H

#include "scalarMeshField.H"
#include "vectorMeshField.H"

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
    #include "symmStencilMeshField.C"
#endif

#endif
