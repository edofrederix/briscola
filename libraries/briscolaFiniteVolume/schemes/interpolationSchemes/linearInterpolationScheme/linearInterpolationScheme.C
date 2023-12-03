#include "linearInterpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearInterpolationScheme<Type,MeshType>::linearInterpolationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearInterpolationScheme<Type,MeshType>::linearInterpolationScheme
(
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
linearInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tInterp
    (
        new meshField<FaceSpace<Type>,MeshType>
        (
            "interp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,MeshType>& Interp = tInterp.ref();

    Interp = Zero;

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllCells(Interp, d, i, j, k)
        for (int f = 0; f < 6; f++)
            Interp(d,i,j,k)[f] =
                fwc(d,i,j,k)[f]*field(d,i,j,k)
              + fwn(d,i,j,k)[f]*field(d,labelVector(i,j,k)+faceOffsets[f]);

    return tInterp;
}

}

}

}
