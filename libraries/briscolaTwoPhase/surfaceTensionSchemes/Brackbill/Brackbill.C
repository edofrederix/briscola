#include "Brackbill.H"
#include "twoPhaseModel.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SigmaModel>
Brackbill<SigmaModel>::Brackbill
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
:
    SigmaModel(tpm, dict)
{}

template<class SigmaModel>
Brackbill<SigmaModel>::Brackbill(const Brackbill& sts)
:
    SigmaModel(sts)
{}

template<class SigmaModel>
Brackbill<SigmaModel>::~Brackbill()
{}

template<class SigmaModel>
void Brackbill<SigmaModel>::correct()
{
    SigmaModel::correct();

    colocatedVectorField::operator=
    (
        this->sigma()*this->kappa()*ex::grad(this->alpha())
      / this->tpm_.template meanRho<colocated>()
    );

    // if (this->fvMsh_.structured())
    //     this->stagForce_ = ...
}

}

}

}
