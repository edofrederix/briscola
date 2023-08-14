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
    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    Sys.singular() = true;

    forAll(field, l)
    forAll(field[l], d)
    {
        const meshDirection<faceScalar,MeshType>& fa =
            field.fvMsh().template
            metrics<MeshType>().faceAreas()[l][d];

        const meshDirection<faceScalar,MeshType>& fd =
            field.fvMsh().template
            metrics<MeshType>().faceDeltas()[l][d];

        meshDirection<stencil,MeshType>& A = Sys.A()[l][d];

        A = Zero;

        const meshDirection<faceScalar,MeshType>& lam = lambda[l][d];

        forAllCells(A, i, j, k)
        {
            A(i,j,k) = lam(i,j,k)*fa(i,j,k)*fd(i,j,k);

            A(i,j,k).center() = - neighborSum(A(i,j,k));
        }
    }

    Sys.b() = Zero;

    return tSys;
}

}

}

}
