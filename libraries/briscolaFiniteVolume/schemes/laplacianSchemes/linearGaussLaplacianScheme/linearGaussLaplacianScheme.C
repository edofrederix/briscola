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
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    laplacianScheme<SType,Type,MeshType>(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
linearGaussLaplacianScheme<SType,Type,MeshType>::linearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    laplacianScheme<SType,Type,MeshType>(dictionary(),fvMsh)
{}

template<class SType, class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussLaplacianScheme<SType,Type,MeshType>::exLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tLap
    (
        new meshField<Type,MeshType>
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,MeshType>& Lap = tLap.ref();

    Lap = Zero;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllFaces(Lap, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = fa(d,ijk)[fd*2]*delta(d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(d,ijk)[fd];

        Lap(d,ijk) += value*(field(d,nei) - field(d,ijk))/cv(d,ijk);
        Lap(d,nei) += value*(field(d,ijk) - field(d,nei))/cv(d,nei);
    }

    return tLap;
}

}

}

}
