#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussLaplacianScheme<Type,MeshType>::linearGaussLaplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    laplacianScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearGaussLaplacianScheme<Type,MeshType>::linearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    laplacianScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussLaplacianScheme<Type,MeshType>::laplacian
(
    const meshField<faceScalar,MeshType>& lambda,
    meshField<Type,MeshType>& field
)
{
    const_cast<meshField<faceScalar,MeshType>&>(lambda).restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    meshField<stencil,MeshType>& A = Sys.A();

    A = Zero;
    Sys.singular() = true;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& fd =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    forAllLevels(A, l, d, i, j, k)
    {
        A(l,d,i,j,k) = lambda(l,d,i,j,k)*fa(l,d,i,j,k)*fd(l,d,i,j,k);
        A(l,d,i,j,k).center() = - neighborSum(A(l,d,i,j,k));
    }

    Sys.b() = Zero;

    // const_cast<meshField<faceScalar,MeshType>&>(lambda).makeShallow();

    return tSys;
}

}

}

}
