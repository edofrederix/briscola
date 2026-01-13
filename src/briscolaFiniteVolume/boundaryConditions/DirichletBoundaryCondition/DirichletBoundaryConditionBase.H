#ifndef DirichletBoundaryConditionBase_H
#define DirichletBoundaryConditionBase_H

#include "boundaryCondition.H"
#include "FastPtrList.H"

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
public:

    // Physical type represented by Type on MeshType
    typedef typename outerProduct<Type, typename MeshType::dimType>::type
        physicalType;

protected:

    // Flattened list of boundary values for each level and direction
    FastPtrList<block<Type>> boundaryValues_;

    // Function to initialize uniform boundary values
    void init(const physicalType&);

public:

    // Construct from components
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and zero boundary value
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const zero&
    );

    // Construct from components and uniform boundary value
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const physicalType&
    );

    // Construct from components and boundary values
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const FastPtrList<block<Type>>&
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
    FastPtrList<block<Type>>& boundaryValues()
    {
        return boundaryValues_;
    }

    // Const access to the boundary values
    const FastPtrList<block<Type>>& boundaryValues() const
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
