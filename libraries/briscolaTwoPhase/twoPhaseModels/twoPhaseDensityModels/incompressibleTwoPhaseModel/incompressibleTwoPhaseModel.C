#include "incompressibleTwoPhaseModel.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(incompressibleTwoPhaseModel, 0);

incompressibleTwoPhaseModel::incompressibleTwoPhaseModel
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    twoPhaseModel(fvMsh, dict),
    rho1_(readScalar(dict.lookup("rho1"))),
    rho2_(readScalar(dict.lookup("rho2")))
{}

incompressibleTwoPhaseModel::incompressibleTwoPhaseModel
(
    const incompressibleTwoPhaseModel& tpm
)
:
    twoPhaseModel(tpm),
    rho1_(tpm.rho1_),
    rho2_(tpm.rho2_)
{}

incompressibleTwoPhaseModel::~incompressibleTwoPhaseModel()
{}

void incompressibleTwoPhaseModel::correctMixture()
{
    rhoc_ = rho2_*alpha_ + rho1_*(1.0-alpha_);

    if (rhosPtr_.valid())
    {
        staggeredScalarField& rhos = rhosPtr_();

        const staggeredScalarField alphas
        (
            stagInterp(alpha_)
        );

        rhos = rho2_*alphas + rho1_*(1.0-alphas);
    }
}

}

}

}
