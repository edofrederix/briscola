#ifndef NeumannBoundaryConditionBase_H
#define NeumannBoundaryConditionBase_H

#include "patchBoundaryCondition.H"
#include "FastPtrList.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class NeumannBoundaryConditionBase
:
    public patchBoundaryCondition<Type,MeshType>
{
public:

    // Physical type represented by Type on MeshType
    typedef typename outerProduct<Type, typename MeshType::dimType>::type
        physicalType;

protected:

    // List of boundary gradients for each direction
    FastPtrList<block<Type>> boundaryGradients_;

    // Function to initialize uniform boundary gradients
    void init(const physicalType&);

public:

    // Construct from components
    NeumannBoundaryConditionBase
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and zero boundary gradient
    NeumannBoundaryConditionBase
    (
        const meshLevel<Type,MeshType>&,
        const boundary&,
        const zero&
    );

    // Construct from components and uniform boundary gradient
    NeumannBoundaryConditionBase
    (
        const meshLevel<Type,MeshType>&,
        const boundary&,
        const physicalType&
    );

    // Construct from components and boundary gradients
    NeumannBoundaryConditionBase
    (
        const meshLevel<Type,MeshType>&,
        const boundary&,
        const FastPtrList<block<Type>>&
    );

    // Copy constructor
    NeumannBoundaryConditionBase
    (
        const NeumannBoundaryConditionBase<Type,MeshType>&
    );

    // Copy constructor resetting level reference
    NeumannBoundaryConditionBase
    (
        const NeumannBoundaryConditionBase<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Access to the boundary gradients
    FastPtrList<block<Type>>& boundaryGradients()
    {
        return boundaryGradients_;
    }

    // Const access to the boundary gradients
    const FastPtrList<block<Type>>& boundaryGradients() const
    {
        return boundaryGradients_;
    }

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
