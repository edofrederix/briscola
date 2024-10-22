#include "harmonicViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
harmonicViscosityMixture<BaseModel>::harmonicViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict)
{
    this->mu_.setRestrictionScheme("harmonicFaceAreaWeighted");
}

template<class BaseModel>
harmonicViscosityMixture<BaseModel>::
harmonicViscosityMixture(const harmonicViscosityMixture& tpm)
:
    BaseModel(tpm)
{
    this->mu_.setRestrictionScheme("harmonicFaceAreaWeighted");
}

template<class BaseModel>
harmonicViscosityMixture<BaseModel>::~harmonicViscosityMixture()
{}

template<class BaseModel>
void harmonicViscosityMixture<BaseModel>::correctMixture()
{
    const tmp<meshField<lowerFaceScalar,typename BaseModel::meshType>> talpha
    (
        this->faceAlpha()
    );

    this->mu_ = 1.0/(talpha()/this->mu2_ + (1.0 - talpha())/this->mu1_);

    BaseModel::correctMixture();
}

}

}

}
