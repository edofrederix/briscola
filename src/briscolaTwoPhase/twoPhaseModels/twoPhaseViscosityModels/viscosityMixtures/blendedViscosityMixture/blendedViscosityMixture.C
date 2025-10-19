#include "blendedViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
blendedViscosityMixture<BaseModel>::blendedViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    C_(dict.lookupOrDefault<scalar>("C", 8))
{}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::
blendedViscosityMixture(const blendedViscosityMixture& tpm)
:
    BaseModel(tpm),
    C_(tpm.C_)
{}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::~blendedViscosityMixture()
{}

template<class BaseModel>
void blendedViscosityMixture<BaseModel>::correctMixture()
{
    const faceField<scalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    forAll(this->mu_, fd)
        this->mu_[fd] =
            (this->mu2_ - this->mu1_)/2.0
          * (tanh(C_*atanh((2.0*alpha[fd] - 1.0)*(1-1e-12))) - 1.0)
          + this->mu2_;

    BaseModel::correctMixture();
}

}

}

}
