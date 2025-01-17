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

void twoPhaseModel::correctMixture()
{
    surfaceTensionSchemePtr_->correct();
}

}

}

}
