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

    // Use base class constructors
    using DirichletBoundaryCondition<Type,MeshType>::
        DirichletBoundaryCondition;

    // Construct from components setting zero value
    noSlipBoundaryCondition
    (
        const meshLevel<Type,MeshType>& level,
        const boundary& b
    )
    :
        DirichletBoundaryCondition<Type,MeshType>
        (
            level,
            b,
            Zero
        )
    {}

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new noSlipBoundaryCondition<Type,MeshType>(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<Type,MeshType>> clone
    (
        const meshLevel<Type,MeshType>& level
    ) const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new noSlipBoundaryCondition<Type,MeshType>(*this, level)
        );
    }
};

}

}

}

#endif
