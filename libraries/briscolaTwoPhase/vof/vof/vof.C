#include "vof.H"

#include "addToRunTimeSelectionTable.H"
#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vof, 0);
defineRunTimeSelectionTable(vof, dictionary);

const scalar vof::threshold = 1e-12;

vof::vof(const dictionary& dict, const fvMesh& fvMsh)
:
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
        true
    )
{
    alpha_.setRestrictionScheme("volumeWeighted");
}

vof::vof(const vof& vf)
:
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
        true
    )
{
    alpha_.setRestrictionScheme("volumeWeighted");
}

vof::~vof()
{}

autoPtr<vof> vof::New
(
    const dictionary& dict,
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

void vof::solve(const staggeredScalarField& U)
{
    solve(coloFaceFlux(U)());
}

}

}

}
