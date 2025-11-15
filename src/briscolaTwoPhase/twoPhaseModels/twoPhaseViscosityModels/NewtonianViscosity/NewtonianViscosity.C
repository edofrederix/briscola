#include "NewtonianViscosity.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class DensityModel>
NewtonianViscosity<DensityModel>::NewtonianViscosity
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    DensityModel(fvMsh, dict),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2")))
{}

template<class DensityModel>
NewtonianViscosity<DensityModel>::
NewtonianViscosity(const NewtonianViscosity& tpm)
:
    DensityModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{}

template<class DensityModel>
void NewtonianViscosity<DensityModel>::correct()
{
    DensityModel::correct();
}

}

}

}
