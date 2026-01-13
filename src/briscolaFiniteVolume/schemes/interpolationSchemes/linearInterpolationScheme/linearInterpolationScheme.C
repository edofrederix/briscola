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
tmp<faceField<Type,MeshType>>
linearInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<faceField<Type,MeshType>> tInterp =
        faceField<Type,MeshType>::New
        (
            "interp("+field.name()+")",
            field.fvMsh()
        );

    faceField<Type,MeshType>& interp = tInterp.ref();

    const faceField<scalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const faceField<scalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllFaces(interp, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        interp[fd](d,ijk) =
            field(d,ijk)*fwc[fd](d,ijk)
          + field(d,nei)*fwn[fd](d,ijk);
    }

    return tInterp;
}

}

}

}
