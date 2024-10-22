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
{
    this->mu_.setRestrictionScheme("blendedViscosityMixture");
}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::
blendedViscosityMixture(const blendedViscosityMixture& tpm)
:
    BaseModel(tpm),
    C_(tpm.C_)
{}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::~blendedViscosityMixture()
{
    this->mu_.setRestrictionScheme("blendedViscosityMixture");
}

template<class BaseModel>
void blendedViscosityMixture<BaseModel>::correctMixture()
{
    const tmp<meshField<lowerFaceScalar,typename BaseModel::meshType>> talpha
    (
        this->faceAlpha()
    );

    this->mu_ =
        (this->mu2_ - this->mu1_)/2.0
      * (tanh(C_*atanh((2.0*talpha() - 1.0)*(1-1e-12))) - 1.0)
      + this->mu2_;

    BaseModel::correctMixture();
}

}

}

}
