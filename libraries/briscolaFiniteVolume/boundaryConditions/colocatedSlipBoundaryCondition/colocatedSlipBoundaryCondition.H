#ifndef colocatedSlipBoundaryCondition_H
#define colocatedSlipBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class colocatedSlipBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
    // Disallow default bitwise assignment
    void operator=(const colocatedSlipBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("slip");

    // Construct from components
    colocatedSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    colocatedSlipBoundaryCondition
    (
        const colocatedSlipBoundaryCondition<Type,MeshType>&
    );

    // Copy constructor resetting field reference
    colocatedSlipBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const colocatedSlipBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new colocatedSlipBoundaryCondition<Type,MeshType>(*this)
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
            new colocatedSlipBoundaryCondition<Type,MeshType>(field,*this)
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

    // Eliminated boundary conditions

        // Eliminated ghosts for all boundaries
        virtual bool eliminated(const bool shifted) const { return true; }

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
    #include "colocatedSlipBoundaryCondition.C"
#endif

#endif
