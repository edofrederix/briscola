#include "twoPhaseMultiVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::twoPhaseMultiVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    BaseModel(fvMsh, dict),
    alphas_(2),
    normalSchemes_(2),
    surfaceTensionSchemes_(2),
    vofs_(2),
    CCL_(2)
{
    forAll(alphas_,i)
    {
        alphas_.set
        (
            i,
            new colocatedScalarField
            (
                "alpha"+Foam::name(i),
                fvMsh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                true,
                false
            )
        );

        alphas_[i].setRestrictionScheme("volumeWeighted");
        alphas_[i] = Zero;
    }

    forAll(normalSchemes_, i)
    {
        normalSchemes_.set
        (
            i,
            normalScheme::New
            (
                fvMsh, dict.subDict("normalScheme"), alphas_[i]
            )
        );
    }

    forAll(surfaceTensionSchemes_,i)
    {
        surfaceTensionSchemes_.set
        (
            i,
            surfaceTensionScheme::New
            (
                fvMsh,
                dict.subDict("surfaceTensionScheme"),
                normalSchemes_[i],
                alphas_[i]
            )
        );
    }

    forAll(vofs_,i)
    {
        vofs_.set
        (
            i,
            vof::New
            (
                fvMsh,
                dict.subDict("vof"),
                normalSchemes_[i],
                alphas_[i]
            )
        );
    }

    forAll(CCL_, i)
    {
        CCL_.set
        (
            i,
            CCL::New
            (
                fvMsh,
                dict.found("CCL")
                ? dict.subDict("CCL")
                : dictionary::null,
                alphas_[i]
            )
        );
    }
}

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::twoPhaseMultiVof
(
    const twoPhaseMultiVof& tpm
)
:
    BaseModel(tpm),
    alphas_(tpm.alphas_),
    normalSchemes_(tpm.normalSchemes_),
    surfaceTensionSchemes_(tpm.surfaceTensionSchemes_),
    vofs_(tpm.vofs_),
    CCL_(tpm.CCL_)
{
    forAll(normalSchemes_,i)
    {
        normalSchemes_[i].correct();
    }

    forAll(surfaceTensionSchemes_,i)
    {
        surfaceTensionSchemes_[i].correct();
    }

    BaseModel::correctMixture();
}

template<class BaseModel>
twoPhaseMultiVof<BaseModel>::~twoPhaseMultiVof()
{}

template<class BaseModel>
tmp<colocatedFaceScalarField>
twoPhaseMultiVof<BaseModel>::surfaceTensionFlux()
{
    tmp<colocatedFaceScalarField> tSurfaceTensionFlux
    (
        new colocatedFaceScalarField("surfaceTensionFlux", this->fvMsh_)
    );

    colocatedFaceScalarField& flux = tSurfaceTensionFlux.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    flux = Zero;
    #endif

    forAll(surfaceTensionSchemes_, i)
    {
        flux += static_cast<colocatedFaceScalarField&>
            (
                surfaceTensionSchemes_[i]
            );
    }

    return tSurfaceTensionFlux;
}

template<class BaseModel>
void twoPhaseMultiVof<BaseModel>::correct()
{
    // Correct multiple marker fields
    forAll(vofs_,i)
    {
        vofs_[i].solve(this->coloFaceFlux()());
    }

    // Compute summed alpha field

    this->alpha_ = Zero;

    forAll(alphas_,i)
    {
        this->alpha_ += alphas_[i];
    }

    // Set bounds for summed alpha field

    this->alpha_.correctBoundaryConditions();

    scalarBlock& alpha = this->alpha_.B();

    forAllBlockLinear(alpha, i)
        alpha(i) =
            Foam::min
            (
                Foam::max
                (
                    round(0.5*alpha(i)/vof::threshold)*2.0*vof::threshold,
                    0.0
                ),
                1.0
            );

    // Correct surface tensions
    forAll(surfaceTensionSchemes_,i)
    {
        surfaceTensionSchemes_[i].correct();
    }

    // Correct interface normals
    forAll(normalSchemes_,i)
    {
        normalSchemes_[i].correct();
    }

    BaseModel::correctMixture();
}

}

}

}
