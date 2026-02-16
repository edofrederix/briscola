#ifndef penalizationDirichletImmersedBoundaryCondition_H
#define penalizationDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class penalizationDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("penalizationDirichlet");

    // Constructor
    penalizationDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    penalizationDirichletImmersedBoundaryCondition
    (
        const penalizationDirichletImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    penalizationDirichletImmersedBoundaryCondition
    (
        const penalizationDirichletImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~penalizationDirichletImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new penalizationDirichletImmersedBoundaryCondition<Type,MeshType>
            (
                *this
            )
        );
    }

    // Construct and return a clone for alternative field
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone
    (
        const meshField<Type,MeshType>& field
    ) const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new penalizationDirichletImmersedBoundaryCondition<Type,MeshType>
            (
                *this,
                field
            )
        );
    }

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#ifdef NoRepository
    #include "penalizationDirichletImmersedBoundaryCondition.C"
#endif

#endif
