#include "curvatureScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(curvatureScheme, 0);
defineRunTimeSelectionTable(curvatureScheme, dictionary);

curvatureScheme::curvatureScheme
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    colocatedScalarField("kappa", fvMsh),
    fvMsh_(fvMsh),
    dict_(dict),
    normal_(normal),
    alpha_(alpha)
{}

curvatureScheme::curvatureScheme(const curvatureScheme& s)
:
    colocatedScalarField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    normal_(s.normal_),
    alpha_(s.alpha_)
{}

curvatureScheme::~curvatureScheme()
{}

autoPtr<curvatureScheme> curvatureScheme::New
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
{
    return curvatureScheme::New(tpm.fvMsh(), dict, tpm.normal(), tpm.alpha());
}

autoPtr<curvatureScheme> curvatureScheme::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
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

    return autoPtr<curvatureScheme>
    (
        cstrIter()(fvMsh, dict, normal, alpha)
    );
}

}

}

}
