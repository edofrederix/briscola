#ifndef emptyBoundaryCondition_H
#define emptyBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyBoundaryCondition
:
    public NeumannBoundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("empty");

    // Construct from components
    emptyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    emptyBoundaryCondition
    (
        const emptyBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    emptyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const emptyBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new emptyBoundaryCondition<Type,MeshType>(*this)
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
            new emptyBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return EMPTYBC;
    }
};

}

}

}

#ifdef NoRepository
    #include "emptyBoundaryCondition.C"
#endif

#endif
