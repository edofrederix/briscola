#ifndef stencilField_H
#define stencilField_H

#include "meshFieldSet.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Stencil is a set of mesh fields with 7 component

template<class Type, class MeshType>
class stencilField
:
    public meshFieldSet<Type,MeshType,7>
{
    // Base type
    typedef meshFieldSet<Type,MeshType,7> base;

public:

    // Use from base
    using base::meshFieldSet;
    using base::operator=;
};

}

}

}

#endif
