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
tmp<meshField<faceScalar,MeshType>> vanLeerLimiterScheme<Type,MeshType>::psi
(
    const meshField<faceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    const meshField<faceScalar,MeshType> r(this->r(phi,field));

    return (r + mag(r))/(1 + mag(r));
}

}

}

}
