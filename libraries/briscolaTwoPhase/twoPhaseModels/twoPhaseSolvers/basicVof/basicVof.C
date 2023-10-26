#include "basicVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
basicVof<BaseModel>::basicVof
(
    const IOdictionary& dict,
    const fvMesh& fvMsh
)
:
    BaseModel(dict, fvMsh),
    vfPtr_(vof::New(dict.subDict("vof"), fvMsh)),
    fluxName_(dict.lookupOrDefault<word>("flux", "phi")),
    velocityName_(dict.lookupOrDefault<word>("velocity", "U"))
{}

template<class BaseModel>
basicVof<BaseModel>::basicVof(const basicVof& tpm)
:
    BaseModel(tpm),
    vfPtr_(tpm.vf().clone()),
    fluxName_(tpm.fluxName_),
    velocityName_(tpm.velocityName_)
{}

template<class BaseModel>
basicVof<BaseModel>::~basicVof()
{}

template<class BaseModel>
void basicVof<BaseModel>::correct()
{
    // Update the volume fraction by solving the Vof advection equation

    if
    (
        this->fvMsh_.db().template
        foundObject<colocatedFaceScalarField>(fluxName_)
    )
    {
        const colocatedFaceScalarField& phi =
            this->fvMsh_.db().template
            lookupObject<colocatedFaceScalarField>(fluxName_);

        vfPtr_->solve(phi);
    }
    else
    {
        const staggeredScalarField& U =
            this->fvMsh_.db().template
            lookupObject<staggeredScalarField>(velocityName_);

        vfPtr_->solve(U);
    }

    // Update mixture

    BaseModel::correctMixture(vfPtr_->alpha());
}

}

}

}
