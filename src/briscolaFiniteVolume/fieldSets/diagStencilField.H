#ifndef diagStencilField_H
#define diagStencilField_H

#include "meshFieldSet.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Diagonal stencil is a set of mesh fields with 1 component

template<class Type, class MeshType>
class diagStencilField
:
    public meshFieldSet<Type,MeshType,1>
{
    // Base type
    typedef meshFieldSet<Type,MeshType,1> base;

public:

    // Use from base
    using base::meshFieldSet;
    using base::operator=;
};

}

}

}

#endif
