#ifndef NeumannBoundaryCondition_H
#define NeumannBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class NeumannBoundaryCondition
:
    public boundaryCondition<Type,MeshType>
{
    // Boundary gradients
    List<Type> boundaryGradients_;

    // Disallow default bitwise assignment
    void operator=(const NeumannBoundaryCondition<Type,MeshType>&);

public:

    // Runtime type information
    TypeName("Neumann");

    // Construct from components
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const boundary&
    );

    // Copy constructor
    NeumannBoundaryCondition(const NeumannBoundaryCondition<Type,MeshType>&);

    // Copy constructor resetting field reference
    NeumannBoundaryCondition
    (
        const meshField<Type,MeshType>&,
        const NeumannBoundaryCondition<Type,MeshType>&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<boundaryCondition<Type,MeshType>>
        (
            new NeumannBoundaryCondition<Type,MeshType>(*this)
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
            new NeumannBoundaryCondition<Type,MeshType>(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Prepare the evaluation of the boundary condition for given level
    virtual void prepare(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label);

    // Eliminated boundary conditions

        // Eliminated ghosts for all boundaries
        virtual bool eliminated(const bool shifted) const { return true; }

        // Return the boundary coefficient for given level and direction
        virtual stencil boundaryCoeff(const label, const label);

        // Return the boundary sources for given level and direction
        virtual tmp<block<Type>> boundarySources(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "NeumannBoundaryCondition.C"
#endif

#endif
