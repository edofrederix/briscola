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
    vfPtr_(vof::New(*this, dict.subDict("vof")))
{}

template<class ViscosityModel>
twoPhaseVof<ViscosityModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    ViscosityModel(tpm),
    vfPtr_(tpm.vf().clone())
{
    this->normalSchemePtr_->correct();
    ViscosityModel::correct();
}

template<class ViscosityModel>
void twoPhaseVof<ViscosityModel>::correct()
{
    // Correct the volume fraction
    vfPtr_->solve(this->coloFaceFlux()());

    // Correct the mesh-specific face volume fraction
    ViscosityModel::correctFaceAlpha();

    // Correct the base model
    ViscosityModel::correct();
}

}

}

}
