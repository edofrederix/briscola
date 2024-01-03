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

    colocatedVectorField::operator=
    (
        this->sigma()*this->kappa()*ex::grad(this->alpha())
      / this->tpm_.template meanRho<colocated>()
    );

    colocatedScalarField meanRho("meanRho",this->tpm_.template meanRho<colocated>()());

    this->surfaceTensionPotential_ =
        ex::interp(this->kappa())
      * ex::faceGrad(this->alpha())
      * ex::interp(this->sigma())
      / ex::interp(meanRho);

    if (this->fvMsh_.structured())
    {
        staggeredScalarField& stagTension = this->stagForcePtr_();
        colocatedLowerFaceScalarField kappaf = ex::interp(this->kappa());
        stagTension = stagInterp(this->sigma())*stagInterp(kappaf)*ex::stagGrad(this->alpha())
                    / this->tpm_.template meanRho<staggered>();
    }
}

}

}

}
