#include "volumeWeightedViscosity.H"
#include "interpolationScheme.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
volumeWeightedViscosity<BaseModel>::volumeWeightedViscosity
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2")))
{
    this->mu_.setRestrictionScheme("faceAreaWeighted");
}

template<class BaseModel>
volumeWeightedViscosity<BaseModel>::volumeWeightedViscosity
(
    const volumeWeightedViscosity& tpm
)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{
    this->mu_.setRestrictionScheme("faceAreaWeighted");
}

template<class BaseModel>
volumeWeightedViscosity<BaseModel>::~volumeWeightedViscosity()
{}

template<class BaseModel>
void volumeWeightedViscosity<BaseModel>::correctMixture()
{
    const meshField<lowerFaceScalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    this->mu_ = alpha*mu2_ + (1.0-alpha)*mu1_;

    BaseModel::correctMixture();
}

}

}

}
