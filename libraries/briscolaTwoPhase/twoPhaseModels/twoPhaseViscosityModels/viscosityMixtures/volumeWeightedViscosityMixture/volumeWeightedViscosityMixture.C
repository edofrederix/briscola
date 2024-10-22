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
    const tmp<meshField<lowerFaceScalar,typename BaseModel::meshType>> talpha
    (
        this->faceAlpha()
    );

    this->mu_ = talpha()*this->mu2_ + (1.0-talpha())*this->mu1_;

    BaseModel::correctMixture();
}

}

}

}
