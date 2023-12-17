#ifndef slipBoundaryCondition_H
#define slipBoundaryCondition_H

#include "boundaryCondition.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class slipBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
public:

    // Runtime type information
    TypeName("slip");

    // Construct from components
    slipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    slipBoundaryCondition
    (
        const slipBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    slipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const slipBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new slipBoundaryCondition<Type,MeshType>(*this)
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
            new slipBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Eliminated boundary condition
    virtual bool eliminated() const { return false; }

    // Evaluate the boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "slipBoundaryCondition.C"
#endif

#endif
