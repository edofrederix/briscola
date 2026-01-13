#include "TwoPhaseModel.H"
#include "exSchemesFaceFlux.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

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

template<>
tmp<colocatedScalarField> TwoPhaseModel<colocated>::meshTypeAlpha() const
{
    return alpha_;
}

template<>
tmp<staggeredScalarField> TwoPhaseModel<staggered>::meshTypeAlpha() const
{
    tmp<staggeredScalarField> tAlpha =
        staggeredScalarField::New("alpha", this->fvMsh_);

    staggeredScalarField& alpha = tAlpha.ref();
    alpha.makeDeep();

    forAllCells(alpha, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        alpha(l,d,i,j,k) =
            0.5*(this->alpha_(l,0,ijk) + this->alpha_(l,0,nei));
    }

    alpha.correctBoundaryConditions();

    return tAlpha;
}

template<>
void TwoPhaseModel<colocated>::correctFaceAlpha()
{
    const colocatedScalarFaceField& fwc =
        this->fvMsh_.template metrics<colocated>().faceWeightsCenter();

    const colocatedScalarFaceField& fwn =
        this->fvMsh_.template metrics<colocated>().faceWeightsNeighbor();

    const colocatedScalarField& alpha = this->alpha_;

    forAllFaces(faceAlpha_, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        faceAlpha_[fd](l,d,ijk) =
            alpha(l,d,ijk)*fwc[fd](l,d,ijk)
          + alpha(l,d,nei)*fwn[fd](l,d,ijk);
    }
}

template<>
void TwoPhaseModel<staggered>::correctFaceAlpha()
{
    const colocatedScalarField& alpha = this->alpha_;

    forAllFaces(faceAlpha_, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        if (d == fd)
        {
            faceAlpha_[fd](l,d,ijk) = alpha(l,0,nei);
        }
        else
        {
            faceAlpha_[fd](l,d,ijk) =
                0.25
              * (
                    alpha(l,0,ijk)
                  + alpha(l,0,lowerNeighbor(ijk,fd))
                  + alpha(l,0,nei)
                  + alpha(l,0,lowerNeighbor(nei,fd))
                );
        }
    }
}

template<class MeshType>
TwoPhaseModel<MeshType>::TwoPhaseModel
(
    const fvMesh& fvMsh,
    const IOdictionary& dict
)
:
    twoPhaseModel(fvMsh, dict),
    faceAlpha_
    (
        "alpha",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    rho_
    (
        "rho",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    rhoMean_(0.0),
    mu_
    (
        "mu",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    )
{
    faceAlpha_ = Zero;
    rho_ = Zero;
    mu_ = Zero;
}

template<class MeshType>
TwoPhaseModel<MeshType>::TwoPhaseModel(const TwoPhaseModel& tpm)
:
    twoPhaseModel(tpm),
    faceAlpha_(tpm.faceAlpha_),
    rho_(tpm.rho_),
    rhoMean_(tpm.rhoMean_),
    mu_(tpm.mu_)
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
tmp<colocatedVectorField> TwoPhaseModel<colocated>::buoyancy() const
{
    tmp<colocatedVectorField> tBuoyancy =
        colocatedVectorField::New("buoyancy", this->fvMsh_);

    colocatedVectorField& buoyancy = tBuoyancy.ref();

    buoyancy = Zero;

    if (Foam::mag(this->g()) > 0.0)
        forAllCells(buoyancy, i, j, k)
            buoyancy(i,j,k) = (rho_(i,j,k) - rhoMean_)*this->g();

    return tBuoyancy;
}

template<>
tmp<staggeredScalarField> TwoPhaseModel<staggered>::buoyancy() const
{
    tmp<staggeredScalarField> tBuoyancy =
        staggeredScalarField::New("buoyancy", this->fvMsh_);

    staggeredScalarField& buoyancy = tBuoyancy.ref();

    buoyancy = Zero;

    if (Foam::mag(this->g()) > 0.0)
    {
        // Project the gravity vector on the staggered base

        const tensor base =
            this->fvMsh_.msh().template cast<rectilinearMesh>().base();

        const vector G = (base & this->g());

        forAllCells(buoyancy, d, i, j, k)
            buoyancy(d,i,j,k) = (rho_(d,i,j,k) - rhoMean_)*G[d];
    }

    return tBuoyancy;
}

// Instantiate
template class TwoPhaseModel<colocated>;
template class TwoPhaseModel<staggered>;

}

}

}
