#include "harmonicViscosityMixture.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
harmonicViscosityMixture<ViscosityModel>::harmonicViscosityMixture
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    ViscosityModel(fvMsh, dict)
{}

template<class ViscosityModel>
harmonicViscosityMixture<ViscosityModel>::
harmonicViscosityMixture(const harmonicViscosityMixture& tpm)
:
    ViscosityModel(tpm)
{}

template<class ViscosityModel>
void harmonicViscosityMixture<ViscosityModel>::correct()
{
    forAllFaces(this->mu_, fd, l, d, i, j, k)
    {
        const scalar alpha = this->faceAlpha_[fd](l,d,i,j,k);

        this->mu_[fd](l,d,i,j,k) =
            1.0/(alpha/this->mu2_ + (1.0 - alpha)/this->mu1_);
    }

    ViscosityModel::correct();
}

}

}

}
