#ifndef faceField_H
#define faceField_H

#include "meshFieldSet.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Face field is a set of mesh fields with 3 components

template<class Type, class MeshType>
class faceField
:
    public meshFieldSet<Type,MeshType,3>
{
    // Base type
    typedef meshFieldSet<Type,MeshType,3> base;

public:

    // Use from base
    using base::meshFieldSet;
    using base::operator=;
};

}

}

}

#endif
