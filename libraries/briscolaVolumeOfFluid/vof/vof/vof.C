#include "vof.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vof, 0);
defineRunTimeSelectionTable(vof, dictionary);

const scalar vof::threshold = 1e-12;

vof::vof(const IOdictionary& dict, const fvMesh& fvMsh)
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
    )
{}

vof::vof(const vof& vf)
:
    regIOobject(vf),
    fvMsh_(vf.fvMsh_),
    dict_(vf.dict_),
    alpha_
    (
        "alpha",
        fvMsh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true,
        false
    )
{}

vof::~vof()
{}

autoPtr<vof> vof::New
(
    const IOdictionary& dict,
    const fvMesh& fvMsh
)
{
    const word vofType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(vofType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Vof scheme type " << vofType
            << ". Valid Vof schemes are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<vof>(cstrIter()(dict, fvMsh));
}

}

}

}
