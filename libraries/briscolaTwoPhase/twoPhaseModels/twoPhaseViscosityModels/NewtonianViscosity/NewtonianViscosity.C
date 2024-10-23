#include "NewtonianViscosity.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
NewtonianViscosity<BaseModel>::NewtonianViscosity
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
NewtonianViscosity<BaseModel>::NewtonianViscosity(const NewtonianViscosity& tpm)
:
    BaseModel(tpm),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_)
{}

template<class BaseModel>
NewtonianViscosity<BaseModel>::~NewtonianViscosity()
{}

template<class BaseModel>
void NewtonianViscosity<BaseModel>::correctMixture()
{
    BaseModel::correctMixture();
}

}

}

}
