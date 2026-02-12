#ifndef dummyBoundaryCondition_H
#define dummyBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class dummyBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("dummy");

    // Construct from components
    dummyBoundaryCondition
    (
        const meshLevel<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    dummyBoundaryCondition(const dummyBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting level reference
    dummyBoundaryCondition
    (
        const dummyBoundaryCondition<Type,MeshType>&,
        const meshLevel<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new dummyBoundaryCondition<Type,MeshType>(*this)
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
            new dummyBoundaryCondition<Type,MeshType>(*this, level)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DUMMYBC;
    }

    // Non-eliminated boundary conditione
    virtual bool eliminated() const { return false; }

    // Evaluate the boundary condition for given direction
    virtual void evaluate(const label);
};

}

}

}

#ifdef NoRepository
    #include "dummyBoundaryCondition.C"
#endif

#endif
