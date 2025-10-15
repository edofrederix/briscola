#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
linearGaussLaplacianScheme<SType,Type,MeshType>::linearGaussLaplacianScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    laplacianScheme<SType,Type,MeshType>(fvMsh, is)
{}

template<class SType, class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussLaplacianScheme<SType,Type,MeshType>::exLaplacian
(
    const meshField<faceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tLap =
        meshField<Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,MeshType>& Lap = tLap.ref();

    const FastPtrList<meshField<scalar,MeshType>>& fa =
        field.fvMsh().template metrics<MeshType>().soa().faceAreas();

    const FastPtrList<meshField<scalar,MeshType>>& delta =
        field.fvMsh().template metrics<MeshType>().soa().faceDeltas();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    #ifdef NO_BLOCK_ZERO_INIT
    Lap = Zero;
    #endif

    forAllFaces(Lap, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value =
            fa[fd](d,ijk)*delta[fd](d,ijk)
          * (field(d,nei) - field(d,ijk))
          * (lambdaPtr ? lambdaPtr->operator()(d,ijk)[fd*2] : 1.0);

        Lap(d,ijk) += value;
        Lap(d,nei) -= value;
    }

    Lap *= icv;

    return tLap;
}

}

}

}
