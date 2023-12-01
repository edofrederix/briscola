#ifndef DirichletBoundaryCondition_H
#define DirichletBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class DirichletBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
    // Boundary values
    List<Type> boundaryValues_;

    // Disallow default bitwise assignment
    void operator=(const DirichletBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("Dirichlet");

    // Construct from components
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    DirichletBoundaryCondition(const DirichletBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const DirichletBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new DirichletBoundaryCondition<Type,MeshType>(*this)
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
            new DirichletBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return DIRICHLETBC;
    }

    // Prepare the evaluation of the boundary condition for given level
    virtual void prepare(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);

    // Constrained boundary conditions

        // This boundary condition constrains values on shifted boundaries
        virtual bool constrained(const bool shifted) const { return shifted; }

        // Return the internal value for given level and direction
        virtual tmp<block<Type>> internalValue(const label, const label);

    // Eliminated boundary conditions

        // Eliminated ghosts for non-shifted boundaries. For shifted boundaries
        // the ghost value is not important.
        virtual bool eliminated(const bool shifted) const { return !shifted; }

        // Return the boundary coefficient for given level and direction
        virtual stencil boundaryCoeff(const label, const label)
        {
            return -pTraits<diagStencil>::one;
        }

        // Return the boundary sources for given level and direction
        virtual tmp<block<Type>> boundarySources(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "DirichletBoundaryCondition.C"
#endif

#endif
