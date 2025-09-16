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

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Construct from components
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Construct from components and boundary values
    DirichletBoundaryConditionBase
    (
        const meshField<Type,MeshType>&,
        const boundary&,
        const List<Type>&
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
