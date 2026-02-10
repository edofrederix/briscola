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
    ViscosityModel(fvMsh, dict)
{
    this->alpha_.init(dict);
}

template<class ViscosityModel>
twoPhaseVof<ViscosityModel>::twoPhaseVof(const twoPhaseVof& tpm)
:
    ViscosityModel(tpm)
{
    this->alpha_.normal().correct();
    this->alpha_.surfaceTension().correct();

    ViscosityModel::correct();
}

template<class ViscosityModel>
void twoPhaseVof<ViscosityModel>::correct()
{
    // Correct the volume fraction
    this->alpha_.vf().solve(this->coloFaceFlux()());

    // Correct the mesh-specific face volume fraction
    ViscosityModel::correctFaceAlpha();

    // Correct the base model
    ViscosityModel::correct();

    // Correct the surface tension
    this->alpha_.surfaceTension().correct();
}

}

}

}
