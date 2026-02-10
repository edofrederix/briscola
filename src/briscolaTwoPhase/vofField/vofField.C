#include "vofField.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vofField, 0);

vofField::vofField(const fvMesh& fvMsh, const word name)
:
    colocatedScalarField
    (
        name == word::null ? "alpha" : name,
        fvMsh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        true,
        true
    ),
    fvMsh_(fvMsh)
{
    colocatedScalarField::setRestrictionScheme("volumeWeighted");

    *this = Zero;
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

void vofField::init(const dictionary& dict)
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
