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

    const faceField<scalar,MeshType> alpha(model.faceAlpha());

    // Set psi values to one for non-interfacial faces

    forAllFaces(psi, fd, d, i, j, k)
        if
        (
            alpha[fd](d,i,j,k) < vof::threshold
         || alpha[fd](d,i,j,k) > 1.0 - vof::threshold
        )
            psi[fd](d,i,j,k) = 1.0;

    return tPsi;
}

}

}

}
