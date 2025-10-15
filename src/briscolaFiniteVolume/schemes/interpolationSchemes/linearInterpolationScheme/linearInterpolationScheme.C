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
    const fvMesh& fvMsh,
    Istream& is
)
:
    interpolationScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
linearInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tInterp =
        meshField<FaceSpace<Type>,MeshType>::New
        (
            "interp("+field.name()+")",
            field.fvMsh()
        );

    meshField<FaceSpace<Type>,MeshType>& Interp = tInterp.ref();

    const faceField<scalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().soa().faceWeightsCenter();

    const faceField<scalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().soa().faceWeightsNeighbor();

    forAllFaces(Interp, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Interp(d,ijk)[fd*2] =
            field(d,ijk)*fwc[fd](d,ijk)
          + field(d,nei)*fwn[fd](d,ijk);

        Interp(d,nei)[fd*2+1] = Interp(d,ijk)[fd*2];
    }

    return tInterp;
}

}

}

}
