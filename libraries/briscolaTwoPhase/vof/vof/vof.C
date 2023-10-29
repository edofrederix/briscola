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

vof::vof
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    colocatedScalarField& alpha
)
:
    regIOobject
    (
        IOobject
        (
            "vof",
            fvMsh.time().timeName(),
            fvMsh.time(),
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

vof::~vof()
{}

autoPtr<vof> vof::New
(
    twoPhaseModel& tpm,
    const dictionary& dict
)
{
    return vof::New(tpm.fvMsh(), dict, tpm.normal(), tpm.alpha());
}

autoPtr<vof> vof::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    colocatedScalarField& alpha
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

void vof::solve(const staggeredScalarField& U)
{
    solve(coloFaceFlux(U)());
}

}

}

}
