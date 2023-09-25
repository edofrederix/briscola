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
    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    forAll(field, l)
    forAll(field[l], d)
    {
        meshDirection<stencil,MeshType>& A = Sys.A()[l][d];

        A = Zero;

        const meshDirection<faceScalar,MeshType>& p = phi[l][d];

        forAllCells(A, i, j, k)
        {
            A(i,j,k) = p(i,j,k);
            forAll(A(i,j,k), cmpt)
            {
                A(i,j,k)[cmpt] = Foam::max(A(i,j,k)[cmpt], 0.0);
            }
            A(i,j,k).center() = neighborSum(A(i,j,k));
        }
    }

    Sys.b() = Zero;

    return tSys;
}

}

}

}
