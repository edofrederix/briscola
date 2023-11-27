#ifndef emptyBoundaryCondition_H
#define emptyBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class emptyBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
private:

    // Disallow default bitwise assignment
    void operator=(const emptyBoundaryCondition<Type,MeshType>&);

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
    emptyBoundaryCondition(const emptyBoundaryCondition<Type,MeshType>&);

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

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);

    // Constrained boundary conditions

        // This boundary condition constrains values on shifted boundaries
        virtual bool constrained(const bool shifted) const { return shifted; }

    // Eliminated boundary conditions

        // Eliminated ghosts for all boundaries
        virtual bool eliminated(const bool shifted) const { return !shifted; }

        // Return the boundary coefficient for given level and direction
        virtual stencil boundaryCoeff(const label l, const label d)
        {
            return pTraits<diagStencil>::one;
        }
};

}

}

}

#ifdef NoRepository
    #include "emptyBoundaryCondition.C"
#endif

#endif
