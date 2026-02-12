#ifndef MittalNeumannImmersedBoundaryCondition_H
#define MittalNeumannImmersedBoundaryCondition_H

#include "MittalImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class MittalNeumannImmersedBoundaryCondition
:
    public MittalImmersedBoundaryCondition<Type,MeshType>
{
    // Boundary gradients
    List<Type> boundaryGradients_;

public:

    // Runtime type information
    TypeName("MittalNeumann");

    // Constructor
    MittalNeumannImmersedBoundaryCondition
    (
        const meshField<Type,MeshType>& field,
        const immersedBoundary<MeshType>& ib
    );

    // Copy constructor
    MittalNeumannImmersedBoundaryCondition
    (
        const MittalNeumannImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    MittalNeumannImmersedBoundaryCondition
    (
        const MittalNeumannImmersedBoundaryCondition&,
        const meshField<Type,MeshType>& field
    );

    // Destructor
    virtual ~MittalNeumannImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<Type,MeshType>> clone() const
    {
        return tmp<immersedBoundaryCondition<Type,MeshType>>
        (
            new MittalNeumannImmersedBoundaryCondition<Type,MeshType>(*this)
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
            new MittalNeumannImmersedBoundaryCondition<Type,MeshType>
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
    #include "MittalNeumannImmersedBoundaryCondition.C"
#endif

#endif
