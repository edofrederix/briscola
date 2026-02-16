#ifndef MittalDirichletImmersedBoundaryCondition_H
#define MittalDirichletImmersedBoundaryCondition_H

#include "MittalImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalDirichletImmersedBoundaryCondition
:
    public MittalImmersedBoundaryCondition<Type,MeshType>
{
    // Boundary values
    List<Type> boundaryValues_;

public:

    // Runtime type information
    TypeName("MittalDirichlet");

    // Constructor
    MittalDirichletImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    MittalDirichletImmersedBoundaryCondition
    (
        const MittalDirichletImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    MittalDirichletImmersedBoundaryCondition
    (
        const MittalDirichletImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~MittalDirichletImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new MittalDirichletImmersedBoundaryCondition<Type,MeshType>(*this)
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
            new MittalDirichletImmersedBoundaryCondition<Type,MeshType>
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
    #include "MittalDirichletImmersedBoundaryCondition.C"
#endif

#endif
