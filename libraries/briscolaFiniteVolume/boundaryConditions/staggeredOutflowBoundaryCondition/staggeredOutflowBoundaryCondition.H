#ifndef staggeredOutflowBoundaryCondition_H
#define staggeredOutflowBoundaryCondition_H

#include "boundaryCondition.H"

#include "block.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class staggeredOutflowBoundaryCondition
:
    public boundaryCondition<scalar,staggered>
{
private:

    // Disallow default bitwise assignment
    void operator=(const staggeredOutflowBoundaryCondition&);

public:

    // Runtime type information
    TypeName("outflow");

    // Construct from components
    staggeredOutflowBoundaryCondition
    (
        const staggeredScalarField&,
        const partPatch&
    );

    // Copy constructor
    staggeredOutflowBoundaryCondition
    (
        const staggeredOutflowBoundaryCondition&
    );

    // Copy constructor resetting field reference
    staggeredOutflowBoundaryCondition
    (
        const staggeredScalarField&,
        const staggeredOutflowBoundaryCondition&
    );

    // Construct and return a clone
    virtual tmp<boundaryCondition<scalar,staggered>> clone() const
    {
        return tmp<boundaryCondition<scalar,staggered>>
        (
            new staggeredOutflowBoundaryCondition(*this)
        );
    }

    // Construct and return a clone setting field reference
    virtual tmp<boundaryCondition<scalar,staggered>> clone
    (
        const staggeredScalarField& field
    ) const
    {
        return tmp<boundaryCondition<scalar,staggered>>
        (
            new staggeredOutflowBoundaryCondition(field,*this)
        );
    }

    // Boundary condition base type
    virtual boundaryConditionBaseType baseType() const
    {
        return NEUMANNBC;
    }

    // Initialize the evaluation of the boundary condition for given level
    virtual void initEvaluate(const label);

    // Evaluate the boundary condition for given level
    virtual void evaluate(const label, const bool = false);
};

}

}

}

#endif
