#ifndef NeumannBoundaryConditionBase_H
#define NeumannBoundaryConditionBase_H

#include "boundaryCondition.H"
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
    public boundaryCondition<Type,MeshType>
{
protected:

    // Flattened list of boundary gradients for each level and direction
    FastPtrList<block<Type>> boundaryGradients_;

    // Function to initialize boundary gradients from a list of uniform
    // gradients for each direction
    void init(const List<Type>&);

public:

    // Construct from components
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and uniform boundary gradient list
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
    );

    // Construct from components and boundary gradients
    NeumannBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const FastPtrList<block<Type>>&
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
