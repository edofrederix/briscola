#include "normalScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(normalScheme, 0);
defineRunTimeSelectionTable(normalScheme, dictionary);

normalScheme::normalScheme(const vof& vf, const dictionary& dict)
:
    vf_(vf),
    dict_(dict)
{}

normalScheme::normalScheme(const normalScheme& s)
:
    vf_(s.vf_),
    dict_(s.dict_)
{}

normalScheme::~normalScheme()
{}

autoPtr<normalScheme> normalScheme::New
(
    const vof& vf,
    const dictionary& dict
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

    return autoPtr<normalScheme>(cstrIter()(vf, dict));
}

}

}

}
