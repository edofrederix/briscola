#include "twoPhaseModels.H"
#include "addToRunTimeSelectionTable.H"

#include "incompressibleTwoPhaseModel.H"
#include "harmonicViscosity.H"
#include "volumeWeightedViscosity.H"
#include "basicVof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Incompressible, harmonic mixing of viscosity and a basic vof solver

makeTwoPhaseModel
(
    basicVof,
    basicVof,
    harmonicViscosity,
    incompressibleTwoPhaseModel
)

// Like the paper of Dodd & Ferrante (2014) with a basic vof solver

makeTwoPhaseModel
(
    basicVofDF,
    basicVof,
    volumeWeightedViscosity,
    incompressibleTwoPhaseModel
)

}

}

}
