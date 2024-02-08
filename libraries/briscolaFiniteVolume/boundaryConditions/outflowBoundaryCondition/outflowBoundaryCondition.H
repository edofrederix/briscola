#ifndef outflowBoundaryCondition_H
#define outflowBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class outflowBoundaryCondition
:
    public NeumannBoundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("outflow");

    // Construct from components
    outflowBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    outflowBoundaryCondition
    (
        const outflowBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    outflowBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const outflowBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new outflowBoundaryCondition<Type,MeshType>(*this)
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
            new outflowBoundaryCondition<Type,MeshType>(field,*this)
        );
    }
};

}

}

}

#ifdef NoRepository
    #include "outflowBoundaryCondition.C"
#endif

#endif
