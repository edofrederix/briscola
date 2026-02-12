#ifndef FadlunDirichletImmersedBoundaryCondition_H
#define FadlunDirichletImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class FadlunDirichletImmersedBoundaryCondition
:
    public immersedBoundaryCondition<Type,MeshType>
{

    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("FadlunDirichlet");

    // Constructor
    FadlunDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    FadlunDirichletImmersedBoundaryCondition
    (
        const FadlunDirichletImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    FadlunDirichletImmersedBoundaryCondition
    (
        const FadlunDirichletImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~FadlunDirichletImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new FadlunDirichletImmersedBoundaryCondition<Type,MeshType>(*this)
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
            new FadlunDirichletImmersedBoundaryCondition<Type,MeshType>
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
    #include "FadlunDirichletImmersedBoundaryCondition.C"
#endif

#endif
