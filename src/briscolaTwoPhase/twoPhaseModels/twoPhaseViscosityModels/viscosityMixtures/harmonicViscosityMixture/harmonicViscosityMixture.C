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
    const auto tAlpha = this->faceAlpha();
    const auto& alpha = tAlpha();

    forAllFaces(this->mu_, fd, d, i, j, k)
        this->mu_[fd](d,i,j,k) =
            1.0
          / (
                alpha[fd](d,i,j,k)/this->mu2_
              + (1.0 - alpha[fd](d,i,j,k))/this->mu1_
            );

    BaseModel::correctMixture();
}

}

}

}
