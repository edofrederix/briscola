#include "volumeWeightedViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::volumeWeightedViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict)
{}

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::volumeWeightedViscosityMixture
(
    const volumeWeightedViscosityMixture& tpm
)
:
    BaseModel(tpm)
{}

template<class BaseModel>
volumeWeightedViscosityMixture<BaseModel>::~volumeWeightedViscosityMixture()
{}

template<class BaseModel>
void volumeWeightedViscosityMixture<BaseModel>::correctMixture()
{
    const auto tAlpha = this->faceAlpha();
    const auto& alpha = tAlpha();

    forAllFaces(this->mu_, fd, d, i, j, k)
        this->mu_[fd](d,i,j,k) =
            alpha[fd](d,i,j,k)*this->mu2_
          + (1.0 - alpha[fd](d,i,j,k))*this->mu1_;

    BaseModel::correctMixture();
}

}

}

}
