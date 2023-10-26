#include "twoPhaseModels.H"
#include "addToRunTimeSelectionTable.H"

#include "incompressibleTwoPhaseModel.H"
#include "harmonicViscosity.H"
#include "basicVof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeTwoPhaseModel
(
    basic,
    basicVof,
    harmonicViscosity,
    incompressibleTwoPhaseModel
)

}

}

}
