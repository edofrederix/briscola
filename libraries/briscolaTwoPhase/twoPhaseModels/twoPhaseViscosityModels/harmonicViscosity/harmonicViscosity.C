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
    const IOdictionary& dict,
    const fvMesh& fvMsh
)
:
    BaseModel(dict, fvMsh),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2")))
{}

template<class BaseModel>
harmonicViscosity<BaseModel>::harmonicViscosity(const harmonicViscosity& tpm)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{}

template<class BaseModel>
harmonicViscosity<BaseModel>::~harmonicViscosity()
{}

template<class BaseModel>
void harmonicViscosity<BaseModel>::correctMixture
(
    const colocatedScalarField& alpha
)
{
    const colocatedFaceScalarField alphaf
    (
        ex::interp(alpha)
    );

    this->muc_ = 1.0/(alphaf/mu2_ + (1.0-alphaf)/mu1_);

    if (this->musPtr_.valid())
    {
        staggeredFaceScalarField& mus = this->musPtr_();

        const staggeredFaceScalarField alphafs
        (
            stagFaceInterp(alpha)
        );

        mus = 1.0/(alphafs/mu2_ + (1.0-alphafs)/mu1_);
    }

    BaseModel::correctMixture(alpha);
}

}

}

}
