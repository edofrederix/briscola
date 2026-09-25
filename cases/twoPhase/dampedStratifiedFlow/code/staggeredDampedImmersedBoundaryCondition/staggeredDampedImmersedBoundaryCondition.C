#include "staggeredDampedImmersedBoundaryCondition.H"
#include "immersedBoundary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(staggeredDampedImmersedBoundaryCondition, 0);

typedef immersedBoundaryCondition<scalar,staggered>
    staggeredScalarImmersedBoundaryCondition;

addToRunTimeSelectionTable
(
    staggeredScalarImmersedBoundaryCondition,
    staggeredDampedImmersedBoundaryCondition,
    dictionary
);

// Constructors

staggeredDampedImmersedBoundaryCondition::
staggeredDampedImmersedBoundaryCondition
(
    const staggeredScalarField& field,
    const immersedBoundary<staggered>& ib
)
:
    immersedBoundaryCondition<scalar,staggered>(field, ib, &ib.mask())
{
    const tensor base =
        this->fvMsh_.msh().template cast<rectilinearMesh>().base();

    d_ =
        (base.x() == vector(unitY) || base.x() == -vector(unitY)) ? 0
      : (base.y() == vector(unitY) || base.y() == -vector(unitY)) ? 1
      : (base.z() == vector(unitY) || base.z() == -vector(unitY)) ? 2
      : -1;

    if (d_ == -1)
        FatalErrorInFunction
            << "The mesh must be aligned with the coordinate system"
            << abort(FatalError);
}

staggeredDampedImmersedBoundaryCondition::
staggeredDampedImmersedBoundaryCondition
(
    const staggeredDampedImmersedBoundaryCondition& ibc
)
:
    immersedBoundaryCondition<scalar,staggered>(ibc),
    d_(ibc.d_)
{}

staggeredDampedImmersedBoundaryCondition::
staggeredDampedImmersedBoundaryCondition
(
    const staggeredDampedImmersedBoundaryCondition& ibc,
    const staggeredScalarField& field
)
:
    immersedBoundaryCondition<scalar,staggered>(ibc, field),
    d_(ibc.d_)
{}

// Destructor

staggeredDampedImmersedBoundaryCondition::
~staggeredDampedImmersedBoundaryCondition()
{}

void staggeredDampedImmersedBoundaryCondition::evaluate
(
    const label l,
    const label d
)
{
    staggeredScalarDirection& x = this->field_[l][d];
    const staggeredScalarDirection& mask = this->forcingMask()[l][d];

    // Set the y-component of the velocity to zero

    if (d == d_)
        forAllCells(x,i,j,k)
            if (mask(i,j,k))
                x(i,j,k) = 0.0;
}

}

}

}
