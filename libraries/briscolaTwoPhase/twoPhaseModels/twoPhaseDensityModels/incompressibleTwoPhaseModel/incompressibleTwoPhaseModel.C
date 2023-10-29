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
    Rho1_(readScalar(dict.lookup("rho1"))),
    Rho2_(readScalar(dict.lookup("rho2")))
{
    this->rho1_ = Rho1_;
    this->rho2_ = Rho2_;
}

incompressibleTwoPhaseModel::incompressibleTwoPhaseModel
(
    const incompressibleTwoPhaseModel& tpm
)
:
    twoPhaseModel(tpm),
    Rho1_(tpm.Rho1_),
    Rho2_(tpm.Rho2_)
{
    this->rho1_ = Rho1_;
    this->rho2_ = Rho2_;
}

incompressibleTwoPhaseModel::~incompressibleTwoPhaseModel()
{}

void incompressibleTwoPhaseModel::correctMixture()
{
    rhoc_ = Rho1_*alpha_ + Rho2_*(1.0-alpha_);

    if (rhosPtr_.valid())
    {
        staggeredScalarField& rhos = rhosPtr_();

        const staggeredScalarField alphas
        (
            stagInterp(alpha_)
        );

        rhos = Rho1_*alphas + Rho2_*(1.0-alphas);
    }
}

}

}

}
