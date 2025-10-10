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
{}

template<class BaseModel>
harmonicViscosityMixture<BaseModel>::
harmonicViscosityMixture(const harmonicViscosityMixture& tpm)
:
    BaseModel(tpm)
{}

template<class BaseModel>
harmonicViscosityMixture<BaseModel>::~harmonicViscosityMixture()
{}

template<class BaseModel>
void harmonicViscosityMixture<BaseModel>::correctMixture()
{
    const meshField<faceScalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    this->mu_ = 1.0/(alpha/this->mu2_ + (1.0 - alpha)/this->mu1_);

    BaseModel::correctMixture();
}

}

}

}
