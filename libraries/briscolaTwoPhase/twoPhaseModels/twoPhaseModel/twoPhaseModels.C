#include "twoPhaseModels.H"
#include "addToRunTimeSelectionTable.H"

#include "TwoPhaseModel.H"
#include "incompressibleTwoPhaseModel.H"
#include "harmonicViscosity.H"
#include "blendedViscosity.H"
#include "volumeWeightedViscosity.H"
#include "twoPhaseVof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Incompressible, blended mixing of viscosity and a basic vof solver

makeTwoPhaseModel
(
    twoPhaseVof,
    twoPhaseVof,
    blendedViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVof,
    twoPhaseVof,
    blendedViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

// Incompressible, harmonic mixing of viscosity and a basic vof solver

makeTwoPhaseModel
(
    twoPhaseVofH,
    twoPhaseVof,
    harmonicViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVofH,
    twoPhaseVof,
    harmonicViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

// Like the paper of Dodd & Ferrante (2014) with a basic vof solver

makeTwoPhaseModel
(
    twoPhaseVofDF,
    twoPhaseVof,
    volumeWeightedViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVofDF,
    twoPhaseVof,
    volumeWeightedViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

}

}

}
