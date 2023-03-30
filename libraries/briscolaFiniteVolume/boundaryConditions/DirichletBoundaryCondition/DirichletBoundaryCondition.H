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
protected:

    // Boundary values for each level and direction
    PtrList<block<Type>> boundaryValues_;

private:

    // Disallow default bitwise assignment
    void operator=(const DirichletBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("Dirichlet");

    // Construct from components
    DirichletBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const partPatch&
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

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label, const bool = false);

    // This boundary condition sets values, so it behaves like a slave patch on
    // shifted boundaries
    virtual bool master() const { return false; }
    virtual bool slave() const { return true; }

    // Eliminated boundary conditions

        // Eliminated ghosts for non-shifted boundaries. For shifted boundaries
        // the values are set directly behaving like a slave boundary.
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
