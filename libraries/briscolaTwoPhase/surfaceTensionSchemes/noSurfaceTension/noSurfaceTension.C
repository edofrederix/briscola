#include "noSurfaceTension.H"
#include "twoPhaseModel.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SigmaModel>
noSurfaceTension<SigmaModel>::noSurfaceTension
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
:
    SigmaModel(tpm, dict)
{
    colocatedLowerFaceScalarField::operator=(Zero);
}

template<class SigmaModel>
noSurfaceTension<SigmaModel>::noSurfaceTension(const noSurfaceTension& sts)
:
    SigmaModel(sts)
{
    colocatedLowerFaceScalarField::operator=(Zero);
}

template<class SigmaModel>
noSurfaceTension<SigmaModel>::~noSurfaceTension()
{}

template<class SigmaModel>
void noSurfaceTension<SigmaModel>::correct()
{
    SigmaModel::correct();
}

}

}

}
