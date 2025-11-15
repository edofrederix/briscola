#include "blendedViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
blendedViscosityMixture<ViscosityModel>::blendedViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    ViscosityModel(fvMsh, dict),
    C_(dict.lookupOrDefault<scalar>("C", 8))
{}

template<class ViscosityModel>
blendedViscosityMixture<ViscosityModel>::
blendedViscosityMixture(const blendedViscosityMixture& tpm)
:
    ViscosityModel(tpm),
    C_(tpm.C_)
{}

template<class ViscosityModel>
void blendedViscosityMixture<ViscosityModel>::correct()
{
    forAllFaces(this->mu_, fd, l, d, i, j, k)
    {
        const scalar alpha = this->faceAlpha_[fd](l,d,i,j,k);

        this->mu_[fd](l,d,i,j,k) =
            (this->mu2_ - this->mu1_)/2.0
          * (
                Foam::tanh
                (
                    C_*Foam::atanh((2.0*alpha - 1.0)*(1-1e-12))
                )
              - 1.0
            )
          + this->mu2_;
    }

    ViscosityModel::correct();
}

}

}

}
