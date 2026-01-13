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
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    SigmaModel(fvMsh, dict, normal, alpha)
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

    colocatedScalarFaceField::operator=
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
