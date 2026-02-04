#include "vofField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vofField, 0);

vofField::vofField
(
    const colocatedScalarField& alpha,
    const dictionary& dict,
    bool initToZero
)
:
    colocatedScalarField
    (
        "alpha",
        alpha,
        false,
        true
    ),
    fvMsh_(alpha.fvMsh()),
    normalSchemePtr_
    (
        normalScheme::New
        (
            fvMsh_,
            dict.subDict("normalScheme"),
            *this
        )
    ),
    surfaceTensionSchemePtr_
    (
        surfaceTensionScheme::New
        (
            fvMsh_,
            dict.subDict("surfaceTensionScheme"),
            normalSchemePtr_(),
            *this
        )
    ),
    vfPtr_
    (
        vof::New
        (
            fvMsh_,
            dict.subDict("vof"),
            normalSchemePtr_(),
            *this
        )
    ),
    tagAlgorithmPtr_
    (
        tagAlgorithm::New
        (
            fvMsh_,
            dict.found("tagAlgorithm")
            ? dict.subDict("tagAlgorithm")
            : dictionary::null,
            *this
        )
    )
{
    colocatedScalarField::setRestrictionScheme("volumeWeighted");

    // Initialize the label field to zero if needed
    if (initToZero)
        static_cast<colocatedScalarField&>(*this) = Zero;
}

vofField::vofField(const fvMesh& fvMsh)
:
    colocatedScalarField
    (
        "alpha",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true
    ),
    fvMsh_(fvMsh),
    normalSchemePtr_(),
    surfaceTensionSchemePtr_(),
    vfPtr_(),
    tagAlgorithmPtr_()
{
    colocatedScalarField::setRestrictionScheme("volumeWeighted");

    static_cast<colocatedScalarField&>(*this) = Zero;
}

vofField::vofField(const vofField& s)
:
    colocatedScalarField(s),
    fvMsh_(s.fvMsh_),
    normalSchemePtr_(s.normalSchemePtr_),
    surfaceTensionSchemePtr_(s.surfaceTensionSchemePtr_),
    vfPtr_(s.vfPtr_),
    tagAlgorithmPtr_(s.tagAlgorithmPtr_)
{}

vofField::~vofField()
{}

void vofField::correctAlpha()
{
    colocatedScalarField::correctBoundaryConditions();

    forAll(*this, l)
    {
        scalarBlock& alpha = (*this)[l].B();

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
    }
}

void vofField::setNormalScheme(const dictionary& dict)
{
    normalSchemePtr_.set
    (
        normalScheme::New
        (
            fvMsh_,
            dict.subDict("normalScheme"),
            *this
        ).ptr()
    );
}

void vofField::setSurfaceTensionScheme(const dictionary& dict)
{
    surfaceTensionSchemePtr_.set
    (
        surfaceTensionScheme::New
        (
            fvMsh_,
            dict.subDict("surfaceTensionScheme"),
            normalSchemePtr_(),
            *this
        ).ptr()
    );
}

void vofField::setVof(const dictionary& dict)
{
    vfPtr_.set
    (
        vof::New
        (
            fvMsh_,
            dict.subDict("vof"),
            normalSchemePtr_(),
            *this
        ).ptr()
    );
}

void vofField::setTagAlgorithm(const dictionary& dict)
{
    tagAlgorithmPtr_.set
    (
        tagAlgorithm::New
        (
            fvMsh_,
            dict.found("tagAlgorithm")
            ? dict.subDict("tagAlgorithm")
            : dictionary::null,
            *this
        ).ptr()
    );
}

}

}

}
