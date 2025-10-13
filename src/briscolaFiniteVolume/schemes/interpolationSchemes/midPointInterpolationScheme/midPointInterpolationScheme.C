#include "midPointInterpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointInterpolationScheme<Type,MeshType>::midPointInterpolationScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    interpolationScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
midPointInterpolationScheme<Type,MeshType>::interp
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

    forAllFaces(Interp, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Interp(d,ijk)[fd*2  ] = 0.5*(field(d,ijk) + field(d,nei));
        Interp(d,nei)[fd*2+1] = Interp(d,ijk)[fd*2];
    }

    return tInterp;
}

}

}

}
