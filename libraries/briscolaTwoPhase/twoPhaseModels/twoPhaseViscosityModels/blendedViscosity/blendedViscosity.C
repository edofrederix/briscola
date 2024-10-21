#include "blendedViscosity.H"
#include "interpolationScheme.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
blendedViscosity<BaseModel>::blendedViscosity
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2"))),
    C_(dict.lookupOrDefault<scalar>("C", 4))
{
    this->mu_.setRestrictionScheme("blendedViscosity");
}

template<class BaseModel>
blendedViscosity<BaseModel>::blendedViscosity(const blendedViscosity& tpm)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_),
    C_(tpm.C_)
{}

template<class BaseModel>
blendedViscosity<BaseModel>::~blendedViscosity()
{
    this->mu_.setRestrictionScheme("blendedViscosity");
}

template<class BaseModel>
void blendedViscosity<BaseModel>::correctMixture()
{
    const meshField<lowerFaceScalar,typename BaseModel::meshType> alpha
    (
        this->faceAlpha()
    );

    this->mu_ =
        (mu2_ - mu1_)/2.0*(tanh(C_*atanh(2.0*alpha - 1.0)) - 1.0) + mu2_;

    BaseModel::correctMixture();
}

}

}

}
