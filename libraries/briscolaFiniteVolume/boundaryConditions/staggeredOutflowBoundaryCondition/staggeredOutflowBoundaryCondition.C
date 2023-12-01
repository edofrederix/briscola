#include "staggeredOutflowBoundaryCondition.H"
#include "addToRunTimeSelectionTable.H"

#include "colocated.H"
#include "staggered.H"

#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(staggeredOutflowBoundaryCondition, 0);

boundaryCondition<scalar,staggered>::
    adddictionaryConstructorToTable<staggeredOutflowBoundaryCondition>
    addstaggeredOutflowBoundaryConditionConstructorToTable_;

staggeredOutflowBoundaryCondition::staggeredOutflowBoundaryCondition
(
    const staggeredScalarField& mshField,
    const boundary& b
)
:
    boundaryCondition<scalar,staggered>(mshField, b)
{}

staggeredOutflowBoundaryCondition::staggeredOutflowBoundaryCondition
(
    const staggeredOutflowBoundaryCondition& bc
)
:
    boundaryCondition<scalar,staggered>(bc.mshField(), bc.mshBoundary())
{}

staggeredOutflowBoundaryCondition::staggeredOutflowBoundaryCondition
(
    const staggeredScalarField& field,
    const staggeredOutflowBoundaryCondition& bc
)
:
    boundaryCondition<scalar,staggered>(field, bc.mshBoundary())
{}

void staggeredOutflowBoundaryCondition::prepare(const label)
{}

void staggeredOutflowBoundaryCondition::evaluate(const label l)
{
    staggeredScalarLevel& field = this->mshField()[l];

    const labelVector bo(this->offset());

    forAll(field, d)
    {
        staggeredScalarDirection& fd = field[d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;

        // Zero-gradient

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = fd(ijk);
        }
    }
}

}

}

}

