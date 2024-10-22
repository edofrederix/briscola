#include "harmonicViscosity.H"
#include "interpolationScheme.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
harmonicViscosity<BaseModel>::harmonicViscosity
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2")))
{
    this->mu_.setRestrictionScheme("harmonicFaceAreaWeighted");
}

template<class BaseModel>
harmonicViscosity<BaseModel>::harmonicViscosity(const harmonicViscosity& tpm)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{
    this->mu_.setRestrictionScheme("harmonicFaceAreaWeighted");
}

template<class BaseModel>
harmonicViscosity<BaseModel>::~harmonicViscosity()
{}

template<class BaseModel>
void harmonicViscosity<BaseModel>::correctMixture()
{
    const meshField<lowerFaceScalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    this->mu_ = 1.0/(alpha/mu2_ + (1.0 - alpha)/mu1_);

    BaseModel::correctMixture();
}

}

}

}
