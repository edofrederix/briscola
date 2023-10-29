#include "surfaceTensionScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(surfaceTensionScheme, 0);
defineRunTimeSelectionTable(surfaceTensionScheme, dictionary);

surfaceTensionScheme::surfaceTensionScheme
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    fvMsh_(fvMsh),
    dict_(dict),
    normal_(normal),
    alpha_(alpha),
    coloForce_
    (
        "surfaceTension",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    sigma_
    (
        "sigma",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    curvatureSchemePtr_
    (
        curvatureScheme::New
        (
            fvMsh,
            dict.subDict("curvatureScheme"),
            normal,
            alpha
        ).ptr()
    )
{
    if (fvMsh.structured())
        stagForcePtr_.reset
        (
            new staggeredScalarField
            (
                "surfaceTension",
                fvMsh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true,
                false
            )
        );
}

surfaceTensionScheme::surfaceTensionScheme(const surfaceTensionScheme& s)
:
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    normal_(s.normal_),
    alpha_(s.alpha_),
    coloForce_(s.coloForce_),
    stagForcePtr_(s.stagForcePtr_, false),
    sigma_(s.sigma_),
    curvatureSchemePtr_(s.curvatureSchemePtr_, false)
{}

surfaceTensionScheme::~surfaceTensionScheme()
{}

autoPtr<surfaceTensionScheme> surfaceTensionScheme::New
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
{
    return surfaceTensionScheme::New
    (
        tpm.fvMsh(),
        dict,
        tpm.normal(),
        tpm.alpha()
    );
}

autoPtr<surfaceTensionScheme> surfaceTensionScheme::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
{
    const word surfaceTensionSchemeType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(surfaceTensionSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown normal scheme type " << surfaceTensionSchemeType
            << ". Valid normal schemes are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<surfaceTensionScheme>
    (
        cstrIter()(fvMsh, dict, normal, alpha)
    );
}

}

}

}
