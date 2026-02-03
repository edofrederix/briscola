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
    dict_(dict),
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

vofField::vofField(const vofField& s)
:
    colocatedScalarField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    normalSchemePtr_(s.normalSchemePtr_),
    surfaceTensionSchemePtr_(s.surfaceTensionSchemePtr_),
    vfPtr_(s.vfPtr_),
    tagAlgorithmPtr_(s.tagAlgorithmPtr_)
{}

vofField::~vofField()
{}

}

}

}
