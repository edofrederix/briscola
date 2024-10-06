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
tmp<meshField<LowerFaceSpace<Type>,MeshType>>
linearInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<LowerFaceSpace<Type>,MeshType>> tInterp
    (
        new meshField<LowerFaceSpace<Type>,MeshType>
        (
            "interp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<LowerFaceSpace<Type>,MeshType>& Interp = tInterp.ref();

    Interp = Zero;

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    forAllFaces(Interp, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Interp(d,ijk)[fd] =
            field(d,ijk)*fwc(d,ijk)[fd*2]
          + field(d,nei)*fwn(d,ijk)[fd*2];
    }

    return tInterp;
}

}

}

}
