#ifndef staggeredDampedImmersedBoundaryCondition_H
#define staggeredDampedImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class staggeredDampedImmersedBoundaryCondition
:
    public immersedBoundaryCondition<scalar,staggered>
{
    // Direction of the y-faces
    label d_;

public:

    // Runtime type information
    TypeName("damped");

    // Constructor
    staggeredDampedImmersedBoundaryCondition
    (
        const staggeredScalarField& field,
        const immersedBoundary<staggered>& ib
    );

    // Copy constructor
    staggeredDampedImmersedBoundaryCondition
    (
        const staggeredDampedImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    staggeredDampedImmersedBoundaryCondition
    (
        const staggeredDampedImmersedBoundaryCondition&,
        const staggeredScalarField& field
    );

    // Destructor
    virtual ~staggeredDampedImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<scalar,staggered>> clone() const
    {
        return tmp<immersedBoundaryCondition<scalar,staggered>>
        (
            new staggeredDampedImmersedBoundaryCondition(*this)
        );
    }

    // Construct and return a clone for alternative field
    virtual tmp<immersedBoundaryCondition<scalar,staggered>> clone
    (
        const staggeredScalarField& field
    ) const
    {
        return tmp<immersedBoundaryCondition<scalar,staggered>>
        (
            new staggeredDampedImmersedBoundaryCondition(*this, field)
        );
    }

    // Evaluate the immersed boundary condition for given level and direction
    virtual void evaluate(const label, const label);
};

}

}

}

#endif
