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
    tmp<meshField<faceScalar,MeshType>> tPsi
    (
        new meshField<faceScalar,MeshType>
        (
            "psi",
            field.fvMsh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false,
            phi.deep() && field.deep()
        )
    );

    meshField<faceScalar,MeshType>& Psi = tPsi.ref();

    Psi = Zero;

    return tPsi;
}

}

}

}
