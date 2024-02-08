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

    rho1_.setRestrictionScheme("volumeWeighted");
    rho2_.setRestrictionScheme("volumeWeighted");
    rhoc_.setRestrictionScheme("volumeWeighted");

    muc_.setRestrictionScheme("harmonicFaceAreaWeighted");

    if (fvMsh_.structured())
    {
        rho1Ptr_->setRestrictionScheme("volumeWeighted");
        rho2Ptr_->setRestrictionScheme("volumeWeighted");
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
    rho1_
    (
        "rho1",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    rho2_
    (
        "rho2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    rhoc_
    (
        "rho",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    muc_
    (
        "mu",
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
    ),
    g_(dict.lookup("g"))
{
    if (fvMsh.structured())
    {
        rho1Ptr_.reset
        (
            new staggeredScalarField
            (
                "rho1",
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true,
                false
            )
        );

        rho2Ptr_.reset
        (
            new staggeredScalarField
            (
                "rho2",
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true,
                false
            )
        );

        rhosPtr_.reset
        (
            new staggeredScalarField
            (
                "rho",
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
            new staggeredLowerFaceScalarField
            (
                "mu",
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
    rho1_(tpm.rho1_),
    rho2_(tpm.rho2_),
    rhoc_(tpm.rhoc_),
    muc_(tpm.muc_),
    rho1Ptr_(tpm.rho1Ptr_, false),
    rho2Ptr_(tpm.rho2Ptr_, false),
    rhosPtr_(tpm.rhosPtr_, false),
    musPtr_(tpm.musPtr_, false),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false),
    g_(tpm.g_)
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
const colocatedLowerFaceScalarField& twoPhaseModel::mu<colocated>() const
{
    return muc_;
}

template<>
const staggeredLowerFaceScalarField& twoPhaseModel::mu<staggered>() const
{
    return musPtr_();
}

template<>
tmp<colocatedScalarField> twoPhaseModel::meanRho<colocated>() const
{
    return (rho1_+rho2_)/2.0;
}

template<>
tmp<staggeredScalarField> twoPhaseModel::meanRho<staggered>() const
{
    return (rho1Ptr_()+rho2Ptr_())/2.0;
}

void twoPhaseModel::correctMixture()
{
    surfaceTensionSchemePtr_->correct();
}

}

}

}
