#include "twoPhaseModel.H"
#include "addToRunTimeSelectionTable.H"

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
    rhoc_.setRestrictionScheme("volumeWeighted");
    muc_.setRestrictionScheme("harmonicFaceAreaWeighted");

    if (fvMsh_.structured())
    {
        rhosPtr_->setRestrictionScheme("volumeWeighted");
        musPtr_->setRestrictionScheme("harmonicFaceAreaWeighted");
    }
}

twoPhaseModel::twoPhaseModel(const IOdictionary& dict, const fvMesh& fvMsh)
:
    regIOobject(dict, true),
    fvMsh_(fvMsh),
    dict_(dict),
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
    rhoc_(tpm.rhoc_),
    muc_(tpm.muc_),
    rhosPtr_(tpm.rhosPtr_, false),
    musPtr_(tpm.musPtr_, false)
{
    setRestrictionSchemes();
}

twoPhaseModel::~twoPhaseModel()
{}

autoPtr<twoPhaseModel> twoPhaseModel::New
(
    const IOdictionary& dict,
    const fvMesh& fvMsh
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

    return autoPtr<twoPhaseModel>(cstrIter()(dict, fvMsh));
}

}

}

}
