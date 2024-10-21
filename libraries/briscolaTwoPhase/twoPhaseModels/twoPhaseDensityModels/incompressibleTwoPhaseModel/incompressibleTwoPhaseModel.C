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

template<class MeshType>
void incompressibleTwoPhaseModel<MeshType>::correctMixture()
{
    const meshField<scalar,MeshType> alpha
    (
        this->meshAlpha()
    );

    this->rho_ = rho2_*alpha + rho1_*(1.0 - alpha);

    TwoPhaseModel<MeshType>::correctMixture();
}

}

}

}
