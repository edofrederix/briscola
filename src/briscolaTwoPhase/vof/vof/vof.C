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

const scalar vof::threshold = 1e-8;

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

void vof::correct()
{
    // Restrict alpha so that derived properties can be computed on all levels

    restrict(alpha_);
    alpha_.correctBoundaryConditions();

    // Apply the alpha correction after the boundary correction and at block
    // level, such that also boundary values are properly set

    forAll(alpha_, l)
    {
        scalarBlock& alpha = alpha_[l].B();

        forAllBlockLinear(alpha, i)
            alpha(i) =
                Foam::min
                (
                    Foam::max
                    (
                        round(0.5*alpha(i)/threshold)*2.0*threshold,
                        0.0
                    ),
                    1.0
                );
    }
}

}

}

}
