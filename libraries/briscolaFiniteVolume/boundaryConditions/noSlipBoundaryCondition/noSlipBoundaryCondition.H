#ifndef noSlipBoundaryCondition_H
#define noSlipBoundaryCondition_H

#include "DirichletBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class noSlipBoundaryCondition
:
    public DirichletBoundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("noSlip");

    // Construct from components
    noSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    noSlipBoundaryCondition
    (
        const noSlipBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    noSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const noSlipBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new noSlipBoundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new noSlipBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DIRICHLETBC;
    }
};

}

}

}

#ifdef NoRepository
    #include "noSlipBoundaryCondition.C"
#endif

#endif
