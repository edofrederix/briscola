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
    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    forAll(field, l)
    forAll(field[l], d)
    {
        meshDirection<stencil,MeshType>& A = Sys.A()[l][d];

        A.initGhosts();

        const meshDirection<faceScalar,MeshType>& p = phi[l][d];
        const meshDirection<Type,MeshType>& f = field[l][d];

        forAllCells(f, i, j, k)
        {
            A(i,j,k) = 0.5*p(i,j,k);

            A(i,j,k).left()   = -A(i,j,k).left();
            A(i,j,k).bottom() = -A(i,j,k).bottom();
            A(i,j,k).aft()    = -A(i,j,k).aft();

            A(i,j,k).center() = neighborSum(A(i,j,k));
        }
    }

    Sys.b() = Zero;

    return tSys;
}

}

}

}
