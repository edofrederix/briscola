#ifndef lubricationBoundaryCondition_H
#define lubricationBoundaryCondition_H

#include "NeumannBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class lubricationBoundaryCondition
:
    public NeumannBoundaryCondition<vector,colocated>
{
    // By default this boundary condition is designed to keep alpha away from
    // the boundary. This behavior can be inverted.
    Switch inverted_;

public:

    // Runtime type information
    TypeName("lubrication");

    // Construct from components
    lubricationBoundaryCondition
    (
        const colocatedVectorLevel&,
        const boundary&
    );

    // Copy constructor
    lubricationBoundaryCondition
    (
        const lubricationBoundaryCondition&
    );

    // Copy constructor resetting level reference
    lubricationBoundaryCondition
    (
        const lubricationBoundaryCondition&,
        const colocatedVectorLevel&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<vector,colocated>> clone() const
    {
        return tmp<boundaryCondition<vector,colocated>>
        (
            new lubricationBoundaryCondition(*this)
        );
    }

    // Construct and return a clone for alternative level
    virtual tmp<boundaryCondition<vector,colocated>> clone
    (
        const colocatedVectorLevel& level
    ) const
    {
        return tmp<boundaryCondition<vector,colocated>>
        (
            new lubricationBoundaryCondition(*this, level)
        );
    }

    // Evaluate the boundary condition for given direction
    virtual void evaluate(const label);
};

}

}

}

#endif
