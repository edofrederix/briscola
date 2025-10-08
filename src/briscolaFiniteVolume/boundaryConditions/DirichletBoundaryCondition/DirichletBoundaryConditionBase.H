#ifndef DirichletBoundaryConditionBase_H
#define DirichletBoundaryConditionBase_H

#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class DirichletBoundaryConditionBase
:
    public boundaryCondition<Type,MeshType>
{
protected:

    // Flattened list of boundary values for each level and direction
    PtrList<block<Type>> boundaryValues_;

    // Function to initialize boundary values from a list of uniform values for
    // each direction
    void init(const List<Type>&);

public:

    // Construct from components
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and uniform boundary value list
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
    );

    // Construct from components and boundary values
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const PtrList<block<Type>>&
    );

    // Copy constructor
    DirichletBoundaryConditionBase
    (
        const DirichletBoundaryConditionBase<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const DirichletBoundaryConditionBase<Type,MeshType>&
    );

    // Access to the boundary values
    PtrList<block<Type>>& boundaryValues()
    {
        return boundaryValues_;
    }

    // Const access to the boundary values
    const PtrList<block<Type>>& boundaryValues() const
    {
        return boundaryValues_;
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DIRICHLETBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "DirichletBoundaryConditionBase.C"
#endif

#endif
