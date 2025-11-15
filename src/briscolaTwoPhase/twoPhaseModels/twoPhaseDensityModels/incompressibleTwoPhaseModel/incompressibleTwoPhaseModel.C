#include "incompressibleTwoPhaseModel.H"
#include "interpolationScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<>
void incompressibleTwoPhaseModel<colocated>::correctRho()
{
    this->rho_ = rho2_*this->alpha_ + rho1_*(1.0 - this->alpha_);
}

template<>
void incompressibleTwoPhaseModel<colocated>::correctMeanRho()
{
    const colocatedScalarDirection& cv =
        this->fvMsh_.template metrics<colocated>().cellVolumes()[0][0];

    const tmp<colocatedScalarField> tMask =
        this->fvMsh_.template metrics<colocated>().fluidMask();

    const colocatedScalarDirection& mask = tMask()[0][0];

    const colocatedScalarDirection& rho = this->rho_[0][0];

    this->rhoMean_ = gSum(mask*rho*cv)/gSum(mask*cv);
}

// Staggered

template<>
void incompressibleTwoPhaseModel<staggered>::correctRho()
{
    const tmp<staggeredScalarField> tAlpha = this->meshTypeAlpha();

    this->rho_ = rho2_*tAlpha() + rho1_*(1.0 - tAlpha());
}

template<>
void incompressibleTwoPhaseModel<staggered>::correctMeanRho()
{
    // Compute the mean mixture mass density on the colocated mesh

    tmp<colocatedScalarDirection> tRho =
        rho2_*this->alpha_[0][0] + rho1_*(1.0 - this->alpha_[0][0]);

    const colocatedScalarDirection& cv =
        this->fvMsh_.template metrics<colocated>().cellVolumes()[0][0];

    const tmp<colocatedScalarField> tMask =
        this->fvMsh_.template metrics<colocated>().fluidMask();

    const colocatedScalarDirection& mask = tMask()[0][0];

    this->rhoMean_ = gSum(mask*tRho()*cv)/gSum(mask*cv);
}

template<class MeshType>
incompressibleTwoPhaseModel<MeshType>::incompressibleTwoPhaseModel
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    TwoPhaseModel<MeshType>(fvMsh, dict),
    rho1_(readScalar(dict.lookup("rho1"))),
    rho2_(readScalar(dict.lookup("rho2")))
{}

template<class MeshType>
incompressibleTwoPhaseModel<MeshType>::incompressibleTwoPhaseModel
(
    const incompressibleTwoPhaseModel& tpm
)
:
    TwoPhaseModel<MeshType>(tpm),
    rho1_(tpm.rho1_),
    rho2_(tpm.rho2_)
{}

template<class MeshType>
void incompressibleTwoPhaseModel<MeshType>::correct()
{
    correctRho();
    correctMeanRho();
    TwoPhaseModel<MeshType>::correct();
}

}

}

}
