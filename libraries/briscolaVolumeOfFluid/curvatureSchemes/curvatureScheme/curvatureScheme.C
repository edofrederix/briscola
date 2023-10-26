#include "curvatureScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(curvatureScheme, 0);
defineRunTimeSelectionTable(curvatureScheme, dictionary);

curvatureScheme::curvatureScheme(const vof& vf, const dictionary& dict)
:
    vf_(vf),
    dict_(dict)
{}

curvatureScheme::curvatureScheme(const curvatureScheme& s)
:
    vf_(s.vf_),
    dict_(s.dict_)
{}

curvatureScheme::~curvatureScheme()
{}

autoPtr<curvatureScheme> curvatureScheme::New
(
    const vof& vf,
    const dictionary& dict
)
{
    const word curvatureSchemeType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(curvatureSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown curvature scheme type " << curvatureSchemeType
            << ". Valid curvature schemes are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<curvatureScheme>(cstrIter()(vf, dict));
}

}

}

}
