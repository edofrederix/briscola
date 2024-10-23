#include "volumeWeightedViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::volumeWeightedViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict)
{
    this->mu_.setRestrictionScheme("faceAreaWeighted");
}

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::volumeWeightedViscosityMixture
(
    const volumeWeightedViscosityMixture& tpm
)
:
    BaseModel(tpm)
{
    this->mu_.setRestrictionScheme("faceAreaWeighted");
}

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::~volumeWeightedViscosityMixture()
{}

template<class BaseModel>
void volumeWeightedViscosityMixture<BaseModel>::correctMixture()
{
    const meshField<lowerFaceScalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    this->mu_ = alpha*this->mu2_ + (1.0-alpha)*this->mu1_;

    BaseModel::correctMixture();
}

}

}

}
