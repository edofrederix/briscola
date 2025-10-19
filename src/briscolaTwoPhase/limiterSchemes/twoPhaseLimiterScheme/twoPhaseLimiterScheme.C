#include "twoPhaseLimiterScheme.H"
#include "TwoPhaseModel.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
twoPhaseLimiterScheme<Type,MeshType>::twoPhaseLimiterScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    limiterScheme<Type,MeshType>(fvMsh, is),
    limiter_
    (
        limiterScheme<Type,MeshType>::New(fvMsh, is)
    )
{}

template<class Type, class MeshType>
tmp<faceField<scalar,MeshType>> twoPhaseLimiterScheme<Type,MeshType>::psi
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    const TwoPhaseModel<MeshType>& model =
        field.fvMsh().db()
       .template lookupObjectRef<twoPhaseModel>("briscolaTwoPhaseDict")
       .template cast<TwoPhaseModel<MeshType>>();

    tmp<faceField<scalar,MeshType>> tPsi = limiter_->psi(phi,field);
    faceField<scalar,MeshType>& psi = tPsi.ref();

    faceField<scalar,MeshType> alpha(model.faceAlpha());

    if (phi[0].deep() && field.deep())
    {
        alpha.setRestrictionScheme("faceAreaWeighted");
        restrict(alpha);
    }

    // Set psi values to zero for non-interfacial faces

    forAllFaces(psi, fd, l, d, i, j, k)
        if
        (
            alpha[fd](l,d,i,j,k) < vof::threshold
         || alpha[fd](l,d,i,j,k) > 1.0 - vof::threshold
        )
            psi[fd](l,d,i,j,k) = 0.0;

    return tPsi;
}

}

}

}
