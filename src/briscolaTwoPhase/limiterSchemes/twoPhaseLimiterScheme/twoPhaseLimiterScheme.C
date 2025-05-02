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
tmp<meshField<faceScalar,MeshType>> twoPhaseLimiterScheme<Type,MeshType>::psi
(
    const meshField<faceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    const TwoPhaseModel<MeshType>& model =
        field.fvMsh().db()
       .template lookupObjectRef<twoPhaseModel>("briscolaTwoPhaseDict")
       .template cast<TwoPhaseModel<MeshType>>();

    const meshField<faceScalar,MeshType> alpha(model.faceAlpha());
    const meshField<faceScalar,MeshType> psi(limiter_->psi(phi,field));

    return
        1.0
      - (1.0 - psi)
      * pos(alpha - vof::threshold)
      * pos(1.0 - vof::threshold - alpha);
}

}

}

}
