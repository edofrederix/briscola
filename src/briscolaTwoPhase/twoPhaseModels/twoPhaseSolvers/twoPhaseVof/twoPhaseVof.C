#include "twoPhaseVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class ViscosityModel>
twoPhaseVof<ViscosityModel>::twoPhaseVof
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    ViscosityModel(fvMsh, dict),
    normalSchemePtr_
    (
        normalScheme::New(*this, dict.subDict("normalScheme")).ptr()
    ),
    surfaceTensionSchemePtr_
    (
        surfaceTensionScheme::New
        (
            fvMsh,
            dict.subDict("surfaceTensionScheme"),
            normalSchemePtr_(),
            this->alpha()
        ).ptr()
    ),
    vfPtr_
    (
        vof::New
        (
            fvMsh,
            dict.subDict("vof"),
            normalSchemePtr_(),
            this->alpha()
        )
    )
{}

template<class ViscosityModel>
twoPhaseVof<ViscosityModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    ViscosityModel(tpm),
    normalSchemePtr_(tpm.normalSchemePtr_, false),
    surfaceTensionSchemePtr_(tpm.surfaceTensionSchemePtr_, false),
    vfPtr_(tpm.vf().clone())
{
    normalSchemePtr_->correct();
    surfaceTensionSchemePtr_->correct();

    ViscosityModel::correct();
}

template<class ViscosityModel>
void twoPhaseVof<ViscosityModel>::correct()
{
    // Correct the volume fraction
    vfPtr_->solve(this->coloFaceFlux()());
    surfaceTensionSchemePtr_->correct();

    ViscosityModel::correct();
}

}

}

}
