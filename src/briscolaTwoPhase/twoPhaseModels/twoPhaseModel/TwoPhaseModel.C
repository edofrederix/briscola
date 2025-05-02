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
tmp<colocatedFaceScalarField> TwoPhaseModel<colocated>::faceAlpha() const
{
    return ex::interp(alpha_);
}

template<>
tmp<staggeredFaceScalarField> TwoPhaseModel<staggered>::faceAlpha() const
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
tmp<colocatedFaceScalarField>
TwoPhaseModel<colocated>::coloFaceFlux() const
{
    return
        fvMsh_.db().template
            lookupObject<colocatedFaceScalarField>("phi");
}

template<>
tmp<colocatedFaceScalarField>
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

template<class MeshType>
tmp<colocatedFaceScalarField> TwoPhaseModel<MeshType>::buoyancy()
{
    setGhA();
    return ghAPtr_()*ex::faceGrad(this->coloRho());
}

template<class MeshType>
tmp<colocatedFaceScalarField> TwoPhaseModel<MeshType>::flux()
{
    tmp<colocatedFaceScalarField> tFlux
    (
        new colocatedFaceScalarField("twoPhaseFlux", this->fvMsh_)
    );

    colocatedFaceScalarField& flux = tFlux.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    flux = Zero;
    #endif

    // And the buoyancy flux when using reduced pressure
    if (this->reduced_)
        flux += this->buoyancy();

    // Add the surface tension flux when we have a surface tension model
    if (this->tension())
        flux +=
            static_cast<colocatedFaceScalarField&>(this->surfaceTension());

    return tFlux;
}

// Instantiate
template class TwoPhaseModel<colocated>;
template class TwoPhaseModel<staggered>;

}

}

}
