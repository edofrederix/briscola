#include "normalScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(normalScheme, 0);
defineRunTimeSelectionTable(normalScheme, dictionary);

normalScheme::normalScheme
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    colocatedVectorField
    (
        "normal",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    alpha_(alpha)
{}

normalScheme::normalScheme(const normalScheme& s)
:
    colocatedVectorField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    alpha_(s.alpha_)
{}

normalScheme::~normalScheme()
{}

autoPtr<normalScheme> normalScheme::New
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
{
    return normalScheme::New(tpm.fvMsh(), dict, tpm.alpha());
}

autoPtr<normalScheme> normalScheme::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
{
    const word normalSchemeType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(normalSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown normal scheme type " << normalSchemeType
            << ". Valid normal schemes are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<normalScheme>(cstrIter()(fvMsh, dict, alpha));
}

}

}

}
