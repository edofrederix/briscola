#include "twoPhaseModel.H"
#include "addToRunTimeSelectionTable.H"
#include "normalScheme.H"
#include "surfaceTensionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPhaseModel, 0);
defineRunTimeSelectionTable(twoPhaseModel, dictionary);

void twoPhaseModel::setRestrictionSchemes()
{
    alpha_.setRestrictionScheme("volumeWeighted");

    rhoc_.setRestrictionScheme("volumeWeighted");
    muc_.setRestrictionScheme("harmonicFaceAreaWeighted");

    if (fvMsh_.structured())
    {
        rhosPtr_->setRestrictionScheme("volumeWeighted");
        musPtr_->setRestrictionScheme("harmonicFaceAreaWeighted");
    }
}

twoPhaseModel::twoPhaseModel(const fvMesh& fvMsh, const IOdictionary& dict)
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
        true,
        false
    ),
    rhoc_
    (
        "rhoc",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    muc_
    (
        "muc",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
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
    )
{
    if (fvMsh.structured())
    {
        rhosPtr_.reset
        (
            new staggeredScalarField
            (
                "rhos",
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true,
                false
            )
        );

        musPtr_.reset
        (
            new staggeredFaceScalarField
            (
                "mus",
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true,
                false
            )
        );
    }

    setRestrictionSchemes();
}

twoPhaseModel::twoPhaseModel(const twoPhaseModel& tpm)
:
    regIOobject(tpm),
    fvMsh_(tpm.fvMsh_),
    dict_(tpm.dict_),
    alpha_(tpm.alpha_),
    rhoc_(tpm.rhoc_),
    muc_(tpm.muc_),
    rhosPtr_(tpm.rhosPtr_, false),
    musPtr_(tpm.musPtr_, false),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false)
{
    setRestrictionSchemes();
}

twoPhaseModel::~twoPhaseModel()
{}

autoPtr<twoPhaseModel> twoPhaseModel::New
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
{
    const word twoPhaseModelType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(twoPhaseModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown two-phase model " << twoPhaseModelType
            << ". Valid two-phase models are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseModel>(cstrIter()(fvMsh, dict));
}

template<>
const colocatedScalarField& twoPhaseModel::rho<colocated>() const
{
    return rhoc_;
}

template<>
const staggeredScalarField& twoPhaseModel::rho<staggered>() const
{
    return rhosPtr_();
}

template<>
const colocatedFaceScalarField& twoPhaseModel::mu<colocated>() const
{
    return muc_;
}

template<>
const staggeredFaceScalarField& twoPhaseModel::mu<staggered>() const
{
    return musPtr_();
}

}

}

}
