#ifndef staggeredSlipBoundaryCondition_H
#define staggeredSlipBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class staggeredSlipBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
    // Disallow default bitwise assignment
    void operator=(const staggeredSlipBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("slip");

    // Construct from components
    staggeredSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
    );

    // Copy constructor
    staggeredSlipBoundaryCondition
    (
        const staggeredSlipBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    staggeredSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const staggeredSlipBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new staggeredSlipBoundaryCondition<Type,MeshType>(*this)
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
            new staggeredSlipBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DIRICHLETBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);

    // Constrained boundary conditions

        // This boundary condition constrains values on shifted boundaries to
        // zero
        virtual bool constrained(const bool shifted) const { return shifted; }

        // Return zero for constrained values
        virtual tmp<block<Type>> internalValue(const label, const label);

    // Eliminated boundary conditions

        // Eliminated ghosts for non-shifted boundaries
        virtual bool eliminated(const bool shifted) const { return !shifted; }

        // Return the boundary coefficient for given level and direction
        virtual stencil boundaryCoeff(const label, const label)
        {
            return pTraits<diagStencil>::one;
        }

        // Return the boundary sources for given level and direction
        virtual tmp<block<Type>> boundarySources(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "staggeredSlipBoundaryCondition.C"
#endif

#endif
