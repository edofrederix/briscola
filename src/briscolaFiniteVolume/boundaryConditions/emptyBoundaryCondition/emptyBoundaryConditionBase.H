#ifndef emptyBoundaryConditionBase_H
#define emptyBoundaryConditionBase_H

#include "boundaryCondition.H"
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
    public boundaryCondition<Type,MeshType>
{
public:

    // Construct from components
    emptyBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    emptyBoundaryConditionBase
    (
        const emptyBoundaryConditionBase<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    emptyBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const emptyBoundaryConditionBase<Type,MeshType>&
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
