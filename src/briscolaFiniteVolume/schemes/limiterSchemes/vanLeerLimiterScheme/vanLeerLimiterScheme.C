#include "vanLeerLimiterScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
vanLeerLimiterScheme<Type,MeshType>::vanLeerLimiterScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    limiterScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<faceField<scalar,MeshType>> vanLeerLimiterScheme<Type,MeshType>::psi
(
    const faceField<scalar,MeshType>& phi,
    const meshField<Type,MeshType>& field,
    const bool deep
)
{
    const faceField<scalar,MeshType> r(this->r(phi,field,deep));
    const faceField<scalar,MeshType> magr(mag(r));

    return (r + magr)/(1.0 + magr);
}

}

}

}
