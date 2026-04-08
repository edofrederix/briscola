#ifndef emptyBoundaryConditionBase_H
#define emptyBoundaryConditionBase_H

#include "patchBoundaryCondition.H"
#include "FastPtrList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyBoundaryConditionBase
:
    public patchBoundaryCondition<Type,MeshType>
{
public:

    // Construct from components
    emptyBoundaryConditionBase
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    emptyBoundaryConditionBase
    (
        const emptyBoundaryConditionBase<Type,MeshType>&
    );

    // Copy constructor resetting level reference
    emptyBoundaryConditionBase
    (
        const emptyBoundaryConditionBase<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return EMPTYBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "emptyBoundaryConditionBase.C"
#endif

#endif
