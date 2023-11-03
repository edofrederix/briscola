#include "upwindDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
upwindDivergenceScheme<Type,MeshType>::upwindDivergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
upwindDivergenceScheme<Type,MeshType>::upwindDivergenceScheme
(
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
upwindDivergenceScheme<Type,MeshType>::div
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

    forAllLevels(A, l, d, i, j, k)
    {
        A(l,d,i,j,k) = phi(l,d,i,j,k);
        forAll(A(i,j,k), cmpt)
        {
            A(i,j,k)[cmpt] = - Foam::min(A(i,j,k)[cmpt], 0.0);
        }
        A(l,d,i,j,k).center() =
            neighborSum(A(l,d,i,j,k)) - neighborSum(phi(l,d,i,j,k));
    }

    Sys.b() = Zero;

    const_cast<meshField<faceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

}

}

}
