#include "twoPhaseModels.H"
#include "addToRunTimeSelectionTable.H"

#include "TwoPhaseModel.H"
#include "incompressibleTwoPhaseModel.H"
#include "harmonicViscosityMixture.H"
#include "blendedViscosityMixture.H"
#include "volumeWeightedViscosityMixture.H"
#include "NewtonianViscosity.H"
#include "twoPhaseVof.H"
#include "twoPhaseMultiVof.H"

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
    blendedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVof,
    twoPhaseVof,
    blendedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

// Incompressible, harmonic mixing of viscosity and a basic vof solver

makeTwoPhaseModel
(
    twoPhaseVofH,
    twoPhaseVof,
    harmonicViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVofH,
    twoPhaseVof,
    harmonicViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

// Like the paper of Dodd & Ferrante (2014) with a basic vof solver

makeTwoPhaseModel
(
    twoPhaseVofDF,
    twoPhaseVof,
    volumeWeightedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseVofDF,
    twoPhaseVof,
    volumeWeightedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

// Incompressible, blended mixing of viscosity and a multiple-marker vof solver

makeTwoPhaseModel
(
    twoPhaseMultiVof,
    twoPhaseMultiVof,
    blendedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    colocated
)

makeTwoPhaseModel
(
    twoPhaseMultiVof,
    twoPhaseMultiVof,
    blendedViscosityMixture,
    NewtonianViscosity,
    incompressibleTwoPhaseModel,
    staggered
)

}

}

}
