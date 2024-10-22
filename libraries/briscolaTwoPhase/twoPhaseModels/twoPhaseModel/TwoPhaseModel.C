#include "TwoPhaseModel.H"
#include "exSchemes.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
TwoPhaseModel<MeshType>::TwoPhaseModel
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    twoPhaseModel(fvMsh, dict),
    rho_
    (
        "rho",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    mu_
    (
        "mu",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    )
{}

template<class MeshType>
TwoPhaseModel<MeshType>::TwoPhaseModel(const TwoPhaseModel& tpm)
:
    twoPhaseModel(tpm),
    rho_(tpm.rho_),
    mu_(tpm.mu_)
{}

template<class MeshType>
TwoPhaseModel<MeshType>::~TwoPhaseModel()
{}

template<class MeshType>
autoPtr<TwoPhaseModel<MeshType>> TwoPhaseModel<MeshType>::New
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
{
    return autoPtr<TwoPhaseModel<MeshType>>
    (
        dynamic_cast<TwoPhaseModel<MeshType>*>
        (
            twoPhaseModel::New<MeshType>(fvMsh, dict).ptr()
        )
    );
}

template<>
tmp<colocatedLowerFaceScalarField> TwoPhaseModel<colocated>::faceAlpha() const
{
    return ex::interp(alpha_);
}

template<>
tmp<staggeredLowerFaceScalarField> TwoPhaseModel<staggered>::faceAlpha() const
{
    return ex::stagFaceInterp(alpha_);
}

template<>
tmp<colocatedScalarField> TwoPhaseModel<colocated>::meshAlpha() const
{
    return alpha_;
}

template<>
tmp<staggeredScalarField> TwoPhaseModel<staggered>::meshAlpha() const
{
    return stagInterp(alpha_);
}

template<>
tmp<colocatedLowerFaceScalarField>
TwoPhaseModel<colocated>::coloFaceFlux() const
{
    return
        fvMsh_.db().template
            lookupObject<colocatedLowerFaceScalarField>("phi");
}

template<>
tmp<colocatedLowerFaceScalarField>
TwoPhaseModel<staggered>::coloFaceFlux() const
{
    return
        ex::coloFaceFlux
        (
            fvMsh_.db().template
                lookupObject<staggeredScalarField>("U")
        );
}

template<>
List<typename colocated::vectorType> TwoPhaseModel<colocated>::gravity() const
{
    return vectorList(1, this->g());
}

template<>
List<typename staggered::vectorType> TwoPhaseModel<staggered>::gravity() const
{
    return list(this->g());
}

}

}

}
