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

template<class MeshType>
scalar TwoPhaseModel<MeshType>::rhoMean() const
{
    const colocatedScalarDirection& cv =
        fvMsh_.template metrics<colocated>().cellVolumes()[0][0];

    return gSum(this->coloRho()()[0][0]*cv)/gSum(cv);
}

template<>
tmp<colocatedVectorField> TwoPhaseModel<colocated>::buoyancy() const
{
    tmp<colocatedVectorField> tBuoyancy
    (
        new colocatedVectorField("buoyancy", this->fvMsh_)
    );

    colocatedVectorField& buoyancy = tBuoyancy.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    buoyancy = Zero;
    #endif

    if (Foam::mag(this->g()) > 0.0)
        buoyancy[0][0] = (this->rho()[0][0] - this->rhoMean())*this->g();

    return tBuoyancy;
}

template<>
tmp<staggeredScalarField> TwoPhaseModel<staggered>::buoyancy() const
{
    tmp<staggeredScalarField> tBuoyancy
    (
        new staggeredScalarField("buoyancy", this->fvMsh_)
    );

    staggeredScalarField& buoyancy = tBuoyancy.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    buoyancy = Zero;
    #endif

    if (Foam::mag(this->g()) > 0.0)
    {
        const scalar rhoMean(this->rhoMean());

        forAll(buoyancy[0], l)
            buoyancy[0][l] = (this->rho()[0][l] - rhoMean)*this->g()[l];
    }

    return tBuoyancy;
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
