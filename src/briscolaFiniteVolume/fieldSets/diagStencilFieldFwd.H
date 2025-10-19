#ifndef diagStencilField_H
#define diagStencilField_H

#include "meshFieldSetFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Diagonal stencil is a set of mesh fields with 1 component

template<class Type, class MeshType>
using diagStencilField = meshFieldSet<Type,MeshType,1>;

}

}

}

#endif
