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
tmp<meshField<LowerFaceSpace<Type>,MeshType>>
midPointInterpolationScheme<Type,MeshType>::interp
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

    forAllFaces(Interp, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(lowerNei(i,j,k,fd));

        Interp(d,ijk)[fd] = 0.5*(field(d,ijk) + field(d,nei));
    }

    return tInterp;
}

}

}

}
