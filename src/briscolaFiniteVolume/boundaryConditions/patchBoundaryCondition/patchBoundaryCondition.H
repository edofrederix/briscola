#ifndef patchBoundaryCondition_H
#define patchBoundaryCondition_H

#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class patchBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
public:

    // Construct from components
    patchBoundaryCondition
    (
        const meshLevel<Type,MeshType>& level,
        const boundary& b
    )
    :
        boundaryCondition<Type,MeshType>(level,b)
    {}

    // Construct from components and zero boundary gradient
    patchBoundaryCondition
    (
        const meshLevel<Type,MeshType>& level,
        const boundary& b,
        const zero&
    )
    :
        boundaryCondition<Type,MeshType>(level,b,Zero)
    {}

    // Copy constructor
    patchBoundaryCondition(const patchBoundaryCondition<Type,MeshType>& bc)
    :
        boundaryCondition<Type,MeshType>(bc)
    {}

    // Copy constructor resetting level reference
    patchBoundaryCondition
    (
        const patchBoundaryCondition<Type,MeshType>& bc,
        const meshLevel<Type,MeshType>& lvl
    )
    :
        boundaryCondition<Type,MeshType>(bc, lvl)
    {}
};

}

}

}

#endif
