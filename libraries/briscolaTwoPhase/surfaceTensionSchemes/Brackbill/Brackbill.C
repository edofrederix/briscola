#include "Brackbill.H"
#include "twoPhaseModel.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

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

    colocatedFaceScalarField::operator=
    (
        this->fvMsh_.template metrics<colocated>().faceAreas()
      * this->kappa().interp()
      * ex::interp(this->sigma_)
      * ex::faceGrad(this->alpha_)
    );
}

}

}

}
