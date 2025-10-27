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
        const meshField<Type,MeshType>& mshField,
        const boundary& b
    )
    :
        DirichletBoundaryCondition<Type,MeshType>
        (
            mshField,
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
};

}

}

}

#endif
