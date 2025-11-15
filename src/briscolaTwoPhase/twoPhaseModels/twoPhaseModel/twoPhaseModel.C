#include "twoPhaseModel.H"
#include "addToRunTimeSelectionTable.H"
#include "normalScheme.H"
#include "surfaceTensionScheme.H"
#include "exSchemesInterpolation.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPhaseModel, 0);

defineRunTimeSelectionTable(twoPhaseModel, colocated);
defineRunTimeSelectionTable(twoPhaseModel, staggered);

twoPhaseModel::twoPhaseModel
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    regIOobject(dict, true),
    fvMsh_(fvMsh),
    dict_(dict),
    alpha_
    (
        "alpha",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    ),
    rho_
    (
        "rho",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    normalSchemePtr_
    (
        normalScheme::New(*this, dict.subDict("normalScheme")).ptr()
    ),
    surfaceTensionSchemePtr_
    (
        surfaceTensionScheme::New
        (
            *this,
            dict.subDict("surfaceTensionScheme")
        ).ptr()
    ),
    g_(dict.lookup("g")),
    tension_(surfaceTensionSchemePtr_->type() != "none")
{
    alpha_.setRestrictionScheme("volumeWeighted");

    alpha_ = Zero;
    rho_ = Zero;
}

twoPhaseModel::twoPhaseModel(const twoPhaseModel& tpm)
:
    regIOobject(tpm),
    fvMsh_(tpm.fvMsh_),
    dict_(tpm.dict_),
    alpha_(tpm.alpha_),
    rho_(tpm.rho_),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false),
    g_(tpm.g_),
    tension_(tpm.tension_)
{
    alpha_.setRestrictionScheme("volumeWeighted");
}

template<>
autoPtr<twoPhaseModel> twoPhaseModel::New<colocated>
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
{
    const word twoPhaseModelType(dict.lookup("type"));

    colocatedConstructorTable::iterator cstrIter =
        colocatedConstructorTablePtr_->find(twoPhaseModelType);

    if (cstrIter == colocatedConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown two-phase model " << twoPhaseModelType
            << ". Valid two-phase models are" << endl
            << colocatedConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseModel>(cstrIter()(fvMsh, dict));
}

template<>
autoPtr<twoPhaseModel> twoPhaseModel::New<staggered>
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
{
    const word twoPhaseModelType(dict.lookup("type"));

    staggeredConstructorTable::iterator cstrIter =
        staggeredConstructorTablePtr_->find(twoPhaseModelType);

    if (cstrIter == staggeredConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown two-phase model " << twoPhaseModelType
            << ". Valid two-phase models are" << endl
            << staggeredConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseModel>(cstrIter()(fvMsh, dict));
}

scalar twoPhaseModel::rhoMean() const
{
    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();

    const tmp<colocatedScalarField> tMask =
        fvMsh_.template metrics<colocated>().fluidMask();

    return gSum(tMask()*rho_*cv)[0]/gSum(tMask()*cv)[0];
}

tmp<colocatedScalarFaceField> twoPhaseModel::flux()
{
    tmp<colocatedScalarFaceField> tFlux
    (
        new colocatedScalarFaceField("twoPhaseFlux", this->fvMsh_)
    );

    colocatedScalarFaceField& flux = tFlux.ref();

    flux = Zero;

    if (this->tension())
        flux +=
            static_cast<colocatedScalarFaceField&>(this->surfaceTension());

    return tFlux;
}

void twoPhaseModel::correct()
{
    surfaceTensionSchemePtr_->correct();
}

}

}

}
