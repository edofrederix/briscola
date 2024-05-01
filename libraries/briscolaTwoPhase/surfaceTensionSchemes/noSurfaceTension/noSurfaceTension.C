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
    colocatedVectorField::operator=(Zero);

    if (this->fvMsh_.structured())
        this->stagForcePtr_() = Zero;
}

template<class SigmaModel>
noSurfaceTension<SigmaModel>::noSurfaceTension(const noSurfaceTension& sts)
:
    SigmaModel(sts)
{
    colocatedVectorField::operator=(Zero);

    this->surfaceTensionPotential_ = Zero;

    if (this->fvMsh_.structured())
        this->stagForcePtr_() = Zero;
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
