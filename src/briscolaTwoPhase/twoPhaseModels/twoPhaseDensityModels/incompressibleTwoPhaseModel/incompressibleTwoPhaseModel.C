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
{
    this->rho_.setRestrictionScheme("volumeWeighted");
}

template<class MeshType>
incompressibleTwoPhaseModel<MeshType>::incompressibleTwoPhaseModel
(
    const incompressibleTwoPhaseModel& tpm
)
:
    TwoPhaseModel<MeshType>(tpm),
    rho1_(tpm.rho1_),
    rho2_(tpm.rho2_)
{
    this->rho_.setRestrictionScheme("volumeWeighted");
}

template<class MeshType>
incompressibleTwoPhaseModel<MeshType>::~incompressibleTwoPhaseModel()
{}

template<>
tmp<colocatedScalarField>
incompressibleTwoPhaseModel<colocated>::coloRho() const
{
    return this->rho_;
}

template<>
tmp<colocatedScalarField>
incompressibleTwoPhaseModel<staggered>::coloRho() const
{
    tmp<colocatedScalarField> tRho
    (
        rho2_*this->alpha_ + rho1_*(1.0 - this->alpha_)
    );

    tRho->correctBoundaryConditions();

    return tRho;
}

template<class MeshType>
void incompressibleTwoPhaseModel<MeshType>::correctMixture()
{
    tmp<meshField<scalar,MeshType>> talpha
    (
        this->meshAlpha()
    );

    this->rho_ = rho2_*talpha() + rho1_*(1.0 - talpha());
    this->rho_.correctBoundaryConditions();

    TwoPhaseModel<MeshType>::correctMixture();
}

}

}

}
