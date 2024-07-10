#include "curvatureScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseModel.H"
#include "vof.H"

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
    colocatedScalarField
    (
        "kappa",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    normal_(normal),
    alpha_(alpha)
{
    // Initialize the curvature to zero

    static_cast<colocatedScalarField&>(*this) = Zero;
}

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

tmp<colocatedLowerFaceScalarField> curvatureScheme::interp() const
{
    const colocatedScalarField& kappa = *this;

    tmp<colocatedLowerFaceScalarField> tInterp
    (
        new colocatedLowerFaceScalarField
        (
            this->name(),
            fvMsh_
        )
    );

    colocatedLowerFaceScalarField& Interp = tInterp.ref();

    Interp = Zero;

    forAllFaces(Interp, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNei(ijk,fd));

        const bool lowerHasInterface =
            alpha_(nei) > vof::threshold && alpha_(nei) < 1 - vof::threshold;

        const bool upperHasInterface =
            alpha_(ijk) > vof::threshold && alpha_(ijk) < 1 - vof::threshold;

        if (lowerHasInterface && upperHasInterface)
        {
            Interp(ijk)[fd] = 0.5*(kappa(ijk) + kappa(nei));
        }
        else if (lowerHasInterface)
        {
            Interp(ijk)[fd] = kappa(nei);
        }
        else if (upperHasInterface)
        {
            Interp(ijk)[fd] = kappa(ijk);
        }
    }

    return tInterp;
}

}

}

}
