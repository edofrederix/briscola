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
    const twoPhaseModel& tpm,
    const dictionary& dict
)
:
    colocatedVectorField
    (
        "surfaceTension",
        tpm.fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    tpm_(tpm),
    fvMsh_(tpm.fvMsh()),
    dict_(dict),
    normal_(tpm.normal()),
    alpha_(tpm.alpha()),
    sigma_
    (
        "sigma",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    sigmaValue_(readScalar(dict.lookup("sigma"))),
    curvatureSchemePtr_
    (
        word(dict.lookup("type")) != "none"
      ? curvatureScheme::New
        (
            fvMsh_,
            dict.subDict("curvatureScheme"),
            normal_,
            alpha_
        ).ptr()
      : nullptr
    ),
    surfaceTensionPotential_
    (
        "surfPot",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    )
{
    if (fvMsh_.structured())
        stagForcePtr_.reset
        (
            new staggeredScalarField
            (
                "surfaceTension",
                fvMsh_,
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
    colocatedVectorField(s),
    tpm_(s.tpm_),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    normal_(s.normal_),
    alpha_(s.alpha_),
    stagForcePtr_(s.stagForcePtr_, false),
    sigma_(s.sigma_),
    sigmaValue_(s.sigmaValue_),
    curvatureSchemePtr_(s.curvatureSchemePtr_, false),
    surfaceTensionPotential_(s.surfaceTensionPotential_)
{}

surfaceTensionScheme::~surfaceTensionScheme()
{}

autoPtr<surfaceTensionScheme> surfaceTensionScheme::New
(
    const twoPhaseModel& tpm,
    const dictionary& dict
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
        cstrIter()(tpm, dict)
    );
}

void surfaceTensionScheme::correct()
{
    if (curvatureSchemePtr_.valid())
        curvatureSchemePtr_->correct();
}

}

}

}
