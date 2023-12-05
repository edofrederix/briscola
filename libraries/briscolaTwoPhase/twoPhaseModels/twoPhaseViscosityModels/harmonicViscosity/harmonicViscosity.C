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
    this->scalarMu1_ = mu1_;
    this->scalarMu2_ = mu2_;
}

template<class BaseModel>
harmonicViscosity<BaseModel>::harmonicViscosity(const harmonicViscosity& tpm)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{
    this->scalarMu1_ = mu1_;
    this->scalarMu2_ = mu2_;
}

template<class BaseModel>
harmonicViscosity<BaseModel>::~harmonicViscosity()
{}

template<class BaseModel>
void harmonicViscosity<BaseModel>::correctMixture()
{
    const colocatedFaceScalarField alphaf
    (
        ex::interp(this->alpha_)
    );

    this->muc_ = 1.0/(alphaf/mu2_ + (1.0-alphaf)/mu1_);

    if (this->musPtr_.valid())
    {
        staggeredFaceScalarField& mus = this->musPtr_();

        const staggeredFaceScalarField alphafs
        (
            stagFaceInterp(this->alpha_)
        );

        mus = 1.0/(alphafs/mu2_ + (1.0-alphafs)/mu1_);
    }

    BaseModel::correctMixture();
}

}

}

}
