#ifndef vectorMeshField_H
#define vectorMeshField_H

#include "vector.H"
#include "scalarMeshField.H"

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
    #include "vectorMeshField.C"
#endif

#endif
