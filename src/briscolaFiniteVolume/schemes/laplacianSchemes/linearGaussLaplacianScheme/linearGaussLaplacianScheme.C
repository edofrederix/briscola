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

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    forAllCells(Lap, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Lap(d,i,j,k) = Zero;
        #endif

        for (label f = 0; f < 6; f++)
            Lap(d,i,j,k) +=
                fa(d,i,j,k)[f]*delta(d,i,j,k)[f]
              * (field(d,neighbor(i,j,k,f)) - field(d,i,j,k))
              * (lambdaPtr ? lambdaPtr->operator()(d,i,j,k)[f] : 1.0);

        Lap(d,i,j,k) *= icv(d,i,j,k);
    }

    return tLap;
}

}

}

}
