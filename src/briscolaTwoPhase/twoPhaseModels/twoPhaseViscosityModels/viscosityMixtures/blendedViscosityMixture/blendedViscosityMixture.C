#include "blendedViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
blendedViscosityMixture<BaseModel>::blendedViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    C_(dict.lookupOrDefault<scalar>("C", 8))
{}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::
blendedViscosityMixture(const blendedViscosityMixture& tpm)
:
    BaseModel(tpm),
    C_(tpm.C_)
{}

template<class BaseModel>
blendedViscosityMixture<BaseModel>::~blendedViscosityMixture()
{}

template<class BaseModel>
void blendedViscosityMixture<BaseModel>::correctMixture()
{
    const auto tAlpha = this->faceAlpha();
    const auto& alpha = tAlpha();

    forAllFaces(this->mu_, fd, d, i, j, k)
        this->mu_[fd](d,i,j,k) =
            (this->mu2_ - this->mu1_)/2.0
          * (
                Foam::tanh
                (
                    C_*Foam::atanh((2.0*alpha[fd](d,i,j,k) - 1.0)*(1-1e-12))
                )
              - 1.0
            )
          + this->mu2_;

    BaseModel::correctMixture();
}

}

}

}
