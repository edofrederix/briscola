#ifndef lubricationImmersedBoundaryCondition_H
#define lubricationImmersedBoundaryCondition_H

#include "immersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class lubricationImmersedBoundaryCondition
:
    public immersedBoundaryCondition<vector,colocated>
{
    // By default this immersed boundary condition is designed to keep alpha
    // away from the boundary. This behavior can be inverted.
    Switch inverted_;

public:

    // Runtime type information
    TypeName("lubrication");

    // Constructor
    lubricationImmersedBoundaryCondition
    (
        const colocatedVectorField& field,
        const immersedBoundary<colocated>& ib
    );

    // Copy constructor
    lubricationImmersedBoundaryCondition
    (
        const lubricationImmersedBoundaryCondition&
    );

    // Copy constructor with alternative field
    lubricationImmersedBoundaryCondition
    (
        const lubricationImmersedBoundaryCondition&,
        const colocatedVectorField& field
    );

    // Destructor
    virtual ~lubricationImmersedBoundaryCondition();

    // Construct and return a clone
    virtual tmp<immersedBoundaryCondition<vector,colocated>> clone() const
    {
        return tmp<immersedBoundaryCondition<vector,colocated>>
        (
            new lubricationImmersedBoundaryCondition(*this)
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
            new lubricationImmersedBoundaryCondition
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

#endif
