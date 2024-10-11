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
    this->muc_.setRestrictionScheme("blendedViscosity");

    if (this->musPtr_.valid())
        this->musPtr_->setRestrictionScheme("blendedViscosity");
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
{}

template<class BaseModel>
void blendedViscosity<BaseModel>::correctMixture()
{
    const colocatedLowerFaceScalarField alphaf
    (
        ex::interp(this->alpha_)
    );

    this->muc_ =
        (mu2_ - mu1_)/2.0*(tanh(C_*atanh(2.0*alphaf - 1.0)) - 1.0) + mu2_;

    if (this->musPtr_.valid())
    {
        const staggeredLowerFaceScalarField alphafs
        (
            stagFaceInterp(this->alpha_)
        );

        this->musPtr_() =
            (mu2_ - mu1_)/2.0*(tanh(C_*atanh(2.0*alphafs - 1.0)) - 1.0) + mu2_;
    }

    BaseModel::correctMixture();
}

}

}

}
