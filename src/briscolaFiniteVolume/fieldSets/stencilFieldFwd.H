#ifndef stencilField_H
#define stencilField_H

#include "meshFieldSetFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Stencil is a set of mesh fields with 7 component

template<class Type, class MeshType>
using stencilField = meshFieldSet<Type,MeshType,7>;

}

}

}

#endif
