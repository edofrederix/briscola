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
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointInterpolationScheme<Type,MeshType>::midPointInterpolationScheme
(
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
midPointInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tInterp
    (
        new meshField<FaceSpace<Type>,MeshType>
        (
            "interpolate("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,MeshType>& Interp = tInterp.ref();

    Interp = Zero;

    forAllDirections(Interp, d, i, j, k)
        Interp(d,i,j,k) =
            0.5
          * FaceSpace<Type>
            (
                field(d,i-1,j,k) + field(d,i,j,k),
                field(d,i+1,j,k) + field(d,i,j,k),
                field(d,i,j-1,k) + field(d,i,j,k),
                field(d,i,j+1,k) + field(d,i,j,k),
                field(d,i,j,k-1) + field(d,i,j,k),
                field(d,i,j,k+1) + field(d,i,j,k)
            );

    return tInterp;
}

}

}

}
