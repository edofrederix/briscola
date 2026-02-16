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
tmp<faceField<Type,MeshType>>
midPointInterpolationScheme<Type,MeshType>::interp
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

    if (sigFpeEnabled())
        interp = Zero;

    forAllFaces(interp, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        interp[fd](d,ijk) = 0.5*(field(d,ijk) + field(d,nei));
    }

    return tInterp;
}

}

}

}
