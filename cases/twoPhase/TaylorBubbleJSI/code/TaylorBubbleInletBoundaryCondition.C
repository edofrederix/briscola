#include "TaylorBubbleInletBoundaryCondition.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

typedef TaylorBubbleInletBoundaryCondition<colocated>
    colocatedTaylorBubbleInletBoundaryCondition;

typedef boundaryCondition<vector,colocated>
    colocatedVectorBoundaryCondition;

defineTypeNameAndDebug
(
    colocatedTaylorBubbleInletBoundaryCondition,
    0
);

addToRunTimeSelectionTable
(
    colocatedVectorBoundaryCondition,
    colocatedTaylorBubbleInletBoundaryCondition,
    dictionary
);

TaylorBubbleInletBoundaryCondition<colocated>::
TaylorBubbleInletBoundaryCondition
(
    const meshLevel<vector,colocated>& level,
    const boundary& b
)
:
    DirichletBoundaryCondition<vector,colocated>(level,b,Zero)
{
    // Lookup some values from this->dict() if needed ...
}

void TaylorBubbleInletBoundaryCondition<colocated>::prepare()
{
    // Hardcode inlet value for now

    const vector Ui(0,0,-0.1);

    block<vector>& bv = this->boundaryValues_[0];
    bv = Ui;
}


// Staggered

typedef TaylorBubbleInletBoundaryCondition<staggered>
    staggeredTaylorBubbleInletBoundaryCondition;

typedef boundaryCondition<scalar,staggered>
    staggeredScalarBoundaryCondition;

defineTypeNameAndDebug
(
    staggeredTaylorBubbleInletBoundaryCondition,
    0
);

addToRunTimeSelectionTable
(
    staggeredScalarBoundaryCondition,
    staggeredTaylorBubbleInletBoundaryCondition,
    dictionary
);

TaylorBubbleInletBoundaryCondition<staggered>::
TaylorBubbleInletBoundaryCondition
(
    const meshLevel<scalar,staggered>& level,
    const boundary& b
)
:
    DirichletBoundaryCondition<scalar,staggered>(level,b,Zero)
{
    // Lookup some values from this->dict() if needed ...
}

void TaylorBubbleInletBoundaryCondition<staggered>::prepare()
{
    const tensor base =
        this->fvMsh_.msh().template cast<rectilinearMesh>().base();

    // Hardcode inlet value for now

    const vector Ui(0,0,-0.1);

    for (int d = 0; d < staggered::numberOfDirections; d++)
    {
        block<scalar>& bv = this->boundaryValues_[d];
        bv = staggered::project(Ui, d, base);
    }
}



}

}

}
