#include "midPointGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointGaussDivergenceScheme<Type,MeshType>::midPointGaussDivergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointGaussDivergenceScheme<Type,MeshType>::midPointGaussDivergenceScheme
(
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
midPointGaussDivergenceScheme<Type,MeshType>::div
(
    const meshField<faceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    const_cast<meshField<faceScalar,MeshType>&>(phi).restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    meshField<stencil,MeshType>& A = Sys.A();
    A = Zero;

    forAllCells(A, l, d, i, j, k)
    {
        A(l,d,i,j,k) = 0.5*phi(l,d,i,j,k);
        A(l,d,i,j,k).center() = neighborSum(A(l,d,i,j,k));
    }

    Sys.b() = Zero;

    const_cast<meshField<faceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

}

}

}
