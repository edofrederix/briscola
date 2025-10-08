#ifndef emptyBoundaryConditionBase_H
#define emptyBoundaryConditionBase_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyBoundaryConditionBase
:
    public NeumannBoundaryCondition<Type,MeshType>
{
public:

    // Use constructors from base
    using NeumannBoundaryCondition<Type,MeshType>::
        NeumannBoundaryCondition;

    // Construct from components setting zero gradient
    emptyBoundaryConditionBase
    (
        const meshField<Type,MeshType>& mshField,
        const boundary& b
    )
    :
        NeumannBoundaryCondition<Type,MeshType>
        (
            mshField,
            b,
            List<Type>(MeshType::numberOfDirections, Zero)
        )
    {}

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return EMPTYBC;
    }
};

}

}

}

#endif
