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
    colocatedScalarFaceField
    (
        "surfaceTension",
        fvMsh
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    normal_(normal),
    alpha_(alpha),
    sigma_
    (
        "sigma",
        fvMsh_
    ),
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
    )
{
    static_cast<colocatedScalarFaceField&>(*this) = Zero;
    sigma_ = Zero;
}

surfaceTensionScheme::surfaceTensionScheme(const surfaceTensionScheme& s)
:
    colocatedScalarFaceField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    normal_(s.normal_),
    alpha_(s.alpha_),
    sigma_(s.sigma_),
    curvatureSchemePtr_(s.curvatureSchemePtr_, false)
{}

surfaceTensionScheme::~surfaceTensionScheme()
{}

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
            << "Unknown surface tension scheme type " << surfaceTensionSchemeType
            << ". Valid surface tension schemes are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<surfaceTensionScheme>
    (
        cstrIter()(fvMsh, dict, normal, alpha)
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
