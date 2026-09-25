#include "colocatedDampedImmersedBoundaryCondition.H"
#include "immersedBoundary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(colocatedDampedImmersedBoundaryCondition, 0);

typedef immersedBoundaryCondition<vector,colocated>
    colocatedVectorImmersedBoundaryCondition;

addToRunTimeSelectionTable
(
    colocatedVectorImmersedBoundaryCondition,
    colocatedDampedImmersedBoundaryCondition,
    dictionary
);

// Constructors

colocatedDampedImmersedBoundaryCondition::
colocatedDampedImmersedBoundaryCondition
(
    const colocatedVectorField& field,
    const immersedBoundary<colocated>& ib
)
:
    immersedBoundaryCondition<vector,colocated>(field, ib, &ib.mask())
{}

colocatedDampedImmersedBoundaryCondition::
colocatedDampedImmersedBoundaryCondition
(
    const colocatedDampedImmersedBoundaryCondition& ibc
)
:
    immersedBoundaryCondition<vector,colocated>(ibc)
{}

colocatedDampedImmersedBoundaryCondition::
colocatedDampedImmersedBoundaryCondition
(
    const colocatedDampedImmersedBoundaryCondition& ibc,
    const colocatedVectorField& field
)
:
    immersedBoundaryCondition<vector,colocated>(ibc, field)
{}

// Destructor

colocatedDampedImmersedBoundaryCondition::
~colocatedDampedImmersedBoundaryCondition()
{}

void colocatedDampedImmersedBoundaryCondition::evaluate
(
    const label l,
    const label d
)
{
    colocatedVectorDirection& x = this->field_[l][d];
    const colocatedScalarDirection& mask = this->forcingMask()[l][d];

    // Set the y-component of the velocity to zero
    forAllCells(x,i,j,k)
        if (mask(i,j,k))
            x(i,j,k).y() = 0.0;
}

}

}

}
