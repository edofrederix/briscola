#include "twoPhaseModel.H"
#include "addToRunTimeSelectionTable.H"
#include "normalScheme.H"
#include "surfaceTensionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPhaseModel, 0);
defineRunTimeSelectionTable(twoPhaseModel, dictionary);

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
    alpha_.setRestrictionScheme("volumeWeighted");
    alpha_ = Zero;
}

twoPhaseModel::twoPhaseModel(const twoPhaseModel& tpm)
:
    regIOobject(tpm),
    fvMsh_(tpm.fvMsh_),
    dict_(tpm.dict_),
    alpha_(tpm.alpha_),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false),
    g_(tpm.g_)
{
    alpha_.setRestrictionScheme("volumeWeighted");
}

twoPhaseModel::~twoPhaseModel()
{}

template<class MeshType>
autoPtr<twoPhaseModel> twoPhaseModel::New
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
{
    const word twoPhaseModelType
    (
        MeshType::typeName
      + word(dict.lookup("type"))
    );

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

void twoPhaseModel::correctMixture()
{
    surfaceTensionSchemePtr_->correct();
}

// Instantiate

template autoPtr<twoPhaseModel> twoPhaseModel::New<colocated>
(
    const fvMesh&,
    const IOdictionary&
);

template autoPtr<twoPhaseModel> twoPhaseModel::New<staggered>
(
    const fvMesh&,
    const IOdictionary&
);

}

}

}
