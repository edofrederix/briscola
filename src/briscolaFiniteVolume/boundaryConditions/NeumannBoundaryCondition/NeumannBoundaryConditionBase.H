#ifndef NeumannBoundaryConditionBase_H
#define NeumannBoundaryConditionBase_H

#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class NeumannBoundaryConditionBase
:
    public boundaryCondition<Type,MeshType>
{
protected:

    // Boundary gradients
    List<Type> boundaryGradients_;

public:

    // Construct from components
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and boundary gradient
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
    );

    // Copy constructor
    NeumannBoundaryConditionBase
    (
        const NeumannBoundaryConditionBase<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const NeumannBoundaryConditionBase<Type,MeshType>&
    );

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "NeumannBoundaryConditionBase.C"
#endif

#endif
