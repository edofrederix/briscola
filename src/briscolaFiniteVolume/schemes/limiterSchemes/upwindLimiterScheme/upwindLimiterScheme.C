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
tmp<faceField<scalar,MeshType>> upwindLimiterScheme<Type,MeshType>::psi
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    tmp<faceField<scalar,MeshType>> tPsi =
        faceField<scalar,MeshType>::New("psi", field.fvMsh());

    faceField<scalar,MeshType>& Psi = tPsi.ref();

    forAll(Psi, fd)
        Psi[fd].make(phi[fd].deep() && field.deep());

    Psi = Zero;

    return tPsi;
}

}

}

}
