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
tmp<colocatedScalarFaceField> TwoPhaseModel<colocated>::faceAlpha() const
{
    return ex::interp(alpha_);
}

template<>
tmp<staggeredScalarFaceField> TwoPhaseModel<staggered>::faceAlpha() const
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
tmp<colocatedScalarFaceField>
TwoPhaseModel<colocated>::coloFaceFlux() const
{
    return
        fvMsh_.db().template
            lookupObject<colocatedScalarFaceField>("phi");
}

template<>
tmp<colocatedScalarFaceField>
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

    const colocatedScalarDirection mask
    (
        fvMsh_.template metrics<colocated>().fluidMask()()[0][0]
    );

    return gSum(mask*this->coloRho()()[0][0]*cv)/gSum(mask*cv);
}

template<>
tmp<colocatedVectorField> TwoPhaseModel<colocated>::buoyancy() const
{
    tmp<colocatedVectorField> tBuoyancy =
        colocatedVectorField::New("buoyancy", this->fvMsh_);

    colocatedVectorDirection& buoyancy = tBuoyancy.ref()[0][0];
    const colocatedScalarDirection& rho = this->rho()[0][0];

    #ifdef NO_BLOCK_ZERO_INIT
    buoyancy = Zero;
    #endif

    if (Foam::mag(this->g()) > 0.0)
        buoyancy = (rho - this->rhoMean())*this->g();

    return tBuoyancy;
}

template<>
tmp<staggeredScalarField> TwoPhaseModel<staggered>::buoyancy() const
{
    tmp<staggeredScalarField> tBuoyancy =
        staggeredScalarField::New("buoyancy", this->fvMsh_);

    staggeredScalarLevel& buoyancy = tBuoyancy.ref()[0];
    const staggeredScalarLevel& rho = this->rho()[0];

    #ifdef NO_BLOCK_ZERO_INIT
    buoyancy = Zero;
    #endif

    // The gravity vector is assumed here to be defined in the staggered
    // coordinate system

    if (Foam::mag(this->g()) > 0.0)
    {
        const scalar rhoMean(this->rhoMean());

        forAll(buoyancy, d)
            buoyancy[d] = (rho[d] - rhoMean)*this->g()[d];
    }

    return tBuoyancy;
}

template<class MeshType>
tmp<colocatedScalarFaceField> TwoPhaseModel<MeshType>::flux()
{
    tmp<colocatedScalarFaceField> tFlux
    (
        new colocatedScalarFaceField("twoPhaseFlux", this->fvMsh_)
    );

    colocatedScalarFaceField& flux = tFlux.ref();

    #ifdef NO_BLOCK_ZERO_INIT
    flux = Zero;
    #endif

    if (this->tension())
        flux +=
            static_cast<colocatedScalarFaceField&>(this->surfaceTension());

    return tFlux;
}

// Instantiate
template class TwoPhaseModel<colocated>;
template class TwoPhaseModel<staggered>;

}

}

}
