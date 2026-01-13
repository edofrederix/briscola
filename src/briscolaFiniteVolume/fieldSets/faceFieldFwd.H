#ifndef faceField_H
#define faceField_H

#include "meshFieldSetFwd.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Face field is a set of mesh fields with 3 components

template<class Type, class MeshType>
using faceField = meshFieldSet<Type,MeshType,3>;

}

}

}

#endif
