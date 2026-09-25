#ifndef colocatedDampedImmersedBoundaryCondition_H
#define colocatedDampedImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class colocatedDampedImmersedBoundaryCondition
:
    public immersedBoundaryCondition<vector,colocated>
{
public:

    // Runtime type information
    TypeName("damped");

    // Constructor
    colocatedDampedImmersedBoundaryCondition
    (
        const colocatedVectorField& field,
        const immersedBoundary<colocated>& ib
    );

    // Copy constructor
    colocatedDampedImmersedBoundaryCondition
    (
        const colocatedDampedImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    colocatedDampedImmersedBoundaryCondition
    (
        const colocatedDampedImmersedBoundaryCondition&,
        const colocatedVectorField& field
    );

    // Destructor
    virtual ~colocatedDampedImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<vector,colocated>> clone() const
    {
        return tmp<immersedBoundaryCondition<vector,colocated>>
        (
            new colocatedDampedImmersedBoundaryCondition(*this)
        );
    }

    // Construct and return a clone for alternative field
    virtual tmp<immersedBoundaryCondition<vector,colocated>> clone
    (
        const colocatedVectorField& field
    ) const
    {
        return tmp<immersedBoundaryCondition<vector,colocated>>
        (
            new colocatedDampedImmersedBoundaryCondition(*this, field)
        );
    }

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#endif
