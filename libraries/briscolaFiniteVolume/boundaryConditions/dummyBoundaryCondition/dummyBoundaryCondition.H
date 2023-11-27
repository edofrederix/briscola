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
    // Disallow default bitwise assignment
    void operator=(const dummyBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("dummy");

    // Construct from components
    dummyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    dummyBoundaryCondition(const dummyBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    dummyBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const dummyBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new dummyBoundaryCondition<Type,MeshType>(*this)
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
            new dummyBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DUMMYBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);
};

}

}

}

#ifdef NoRepository
    #include "dummyBoundaryCondition.C"
#endif

#endif
