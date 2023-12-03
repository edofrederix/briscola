#include "linearGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussDivergenceScheme<Type,MeshType>::linearGaussDivergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearGaussDivergenceScheme<Type,MeshType>::linearGaussDivergenceScheme
(
    const fvMesh& fvMsh
)
:
    divergenceScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussDivergenceScheme<Type,MeshType>::div
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

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllCells(A, l, d, i, j, k)
    {
        for (int f = 0; f < 6; f++)
        {
            A(l,d,i,j,k)[0]  += phi(l,d,i,j,k)[f] * fwc(l,d,i,j,k)[f];
            A(l,d,i,j,k)[f+1] = phi(l,d,i,j,k)[f] * fwn(l,d,i,j,k)[f];
        }
    }

    Sys.b() = Zero;

    const_cast<meshField<faceScalar,MeshType>&>(phi).makeShallow();

    return tSys;
}

}

}

}
