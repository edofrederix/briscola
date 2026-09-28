#include "lubricationImmersedBoundaryCondition.H"
#include "immersedBoundary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(lubricationImmersedBoundaryCondition, 0);

typedef immersedBoundaryCondition<vector,colocated>
    colocatedVectorImmersedBoundaryCondition;

addToRunTimeSelectionTable
(
    colocatedVectorImmersedBoundaryCondition,
    lubricationImmersedBoundaryCondition,
    dictionary
);

// Constructors

lubricationImmersedBoundaryCondition::lubricationImmersedBoundaryCondition
(
    const colocatedVectorField& field,
    const immersedBoundary<colocated>& ib
)
:
    immersedBoundaryCondition<vector,colocated>(field, ib, &ib.wallAdjMask()),
    inverted_(this->dict_.lookupOrDefault<Switch>("inverted", false))
{}

lubricationImmersedBoundaryCondition::lubricationImmersedBoundaryCondition
(
    const lubricationImmersedBoundaryCondition& ibc
)
:
    immersedBoundaryCondition<vector,colocated>(ibc),
    inverted_(ibc.inverted_)
{}

lubricationImmersedBoundaryCondition::lubricationImmersedBoundaryCondition
(
    const lubricationImmersedBoundaryCondition& ibc,
    const colocatedVectorField& field
)
:
    immersedBoundaryCondition<vector,colocated>(ibc, field),
    inverted_(ibc.inverted_)
{}

// Destructor

lubricationImmersedBoundaryCondition::~lubricationImmersedBoundaryCondition()
{}

void lubricationImmersedBoundaryCondition::evaluate
(
    const label l,
    const label d
)
{
    colocatedVectorDirection& x = this->field_[l][d];

    const colocatedScalarDirection& mask = this->forcingMask()[l][d];
    const colocatedVectorDirection& normal = this->ib().wallNormalAdj()[l][d];

    // By definition, the interface normal points into regions of high alpha.
    // The immersed boundary condition is designed to keep alpha away from the
    // boundary. So, the interface normal should be set as the wall normal,
    // which always points into the fluid.

    forAllCells(x,i,j,k)
        if (mask(i,j,k))
            x(i,j,k) = normal(i,j,k);
}

}

}

}
