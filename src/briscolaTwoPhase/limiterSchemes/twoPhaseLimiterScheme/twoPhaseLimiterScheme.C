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
    const meshField<Type,MeshType>& field,
    const bool deep
)
{
    const TwoPhaseModel<MeshType>& model =
        field.fvMsh().db()
       .template lookupObjectRef<twoPhaseModel>("briscolaTwoPhaseDict")
       .template cast<TwoPhaseModel<MeshType>>();

    tmp<faceField<scalar,MeshType>> tPsi = limiter_->psi(phi,field,deep);
    faceField<scalar,MeshType>& psi = tPsi.ref();

    const faceField<scalar,MeshType>& faceAlpha = model.faceAlpha();

    // Set psi values to one for non-interfacial faces

    forAllFaces(psi, fd, l, d, i, j, k)
    {
        const scalar alpha = faceAlpha[fd](l,d,i,j,k);

        if (alpha < vof::threshold || alpha > 1.0 - vof::threshold)
        {
            psi[fd](l,d,i,j,k) = 1.0;
        }
    }

    return tPsi;
}

}

}

}
