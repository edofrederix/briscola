#include "upwindLimiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
upwindLimiterScheme<Type,MeshType>::upwindLimiterScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    limiterScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<meshField<faceScalar,MeshType>> upwindLimiterScheme<Type,MeshType>::psi
(
    const meshField<faceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<faceScalar,MeshType>> tPsi =
        meshField<faceScalar,MeshType>::New("psi", field.fvMsh());

    meshField<faceScalar,MeshType>& Psi = tPsi.ref();

    Psi.make(phi.deep() && field.deep());

    Psi = Zero;

    return tPsi;
}

}

}

}
