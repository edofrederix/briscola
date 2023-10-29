#include "twoPhaseVof.H"

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
twoPhaseVof<BaseModel>::twoPhaseVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    vfPtr_(vof::New(*this, dict.subDict("vof"))),
    fluxName_(dict.lookupOrDefault<word>("flux", "phi")),
    velocityName_(dict.lookupOrDefault<word>("velocity", "U"))
{}

template<class BaseModel>
twoPhaseVof<BaseModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    BaseModel(tpm),
    vfPtr_(tpm.vf().clone()),
    fluxName_(tpm.fluxName_),
    velocityName_(tpm.velocityName_)
{}

template<class BaseModel>
twoPhaseVof<BaseModel>::~twoPhaseVof()
{}

template<class BaseModel>
void twoPhaseVof<BaseModel>::correct()
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

    // Correct the mixture

    BaseModel::correctMixture();
}

}

}

}
