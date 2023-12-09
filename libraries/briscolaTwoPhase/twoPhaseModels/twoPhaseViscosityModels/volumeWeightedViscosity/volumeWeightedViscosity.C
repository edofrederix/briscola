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
{}

template<class BaseModel>
volumeWeightedViscosity<BaseModel>::volumeWeightedViscosity
(
    const volumeWeightedViscosity& tpm
)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{}

template<class BaseModel>
volumeWeightedViscosity<BaseModel>::~volumeWeightedViscosity()
{}

template<class BaseModel>
void volumeWeightedViscosity<BaseModel>::correctMixture()
{
    const colocatedLowerFaceScalarField alphaf
    (
        ex::interp(this->alpha_)
    );

    this->muc_ = alphaf*mu2_ + (1.0-alphaf)*mu1_;

    if (this->musPtr_.valid())
    {
        staggeredLowerFaceScalarField& mus = this->musPtr_();

        const staggeredLowerFaceScalarField alphafs
        (
            stagFaceInterp(this->alpha_)
        );

        mus = alphafs*mu2_ + (1.0-alphafs)*mu1_;
    }

    BaseModel::correctMixture();
}

}

}

}
