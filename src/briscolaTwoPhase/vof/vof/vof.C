#include "vof.H"

#include "addToRunTimeSelectionTable.H"
#include "faceFluxScheme.H"

#include "vofField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vof, 0);
defineRunTimeSelectionTable(vof, dictionary);

const scalar vof::threshold = 1e-8;

vof::vof
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    vofField& alpha
)
:
    regIOobject
    (
        IOobject
        (
            "vof",
            fvMsh.time().name(),
            fvMsh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    normal_(normal),
    alpha_(alpha)
{}

vof::vof(const vof& vf)
:
    regIOobject(vf),
    fvMsh_(vf.fvMsh_),
    dict_(vf.dict_),
    normal_(vf.normal_),
    alpha_(vf.alpha_)
{}

autoPtr<vof> vof::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    vofField& alpha
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

    return autoPtr<vof>(cstrIter()(fvMsh, dict, normal, alpha));
}

void vof::correct()
{
    // Restrict alpha so that derived properties can be computed on all levels
    restrict(alpha_);

    // Correct boundary conditions and set value bounds
    alpha_.correctAlpha();
}

}

}

}
