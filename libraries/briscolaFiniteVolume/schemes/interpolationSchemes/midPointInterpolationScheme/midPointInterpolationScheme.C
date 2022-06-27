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

    forAll(field, l)
    forAll(field[l], d)
    {
        meshDirection<FaceSpace<Type>,MeshType>& I = Interp[l][d];
        const meshDirection<Type,MeshType>& f = field[l][d];

        I.initGhosts();

        forAllCells(f, i, j, k)
        {
            I(i,j,k) =
                0.5
              * FaceSpace<Type>
                (
                    f(i-1,j,k) + f(i,j,k),
                    f(i+1,j,k) + f(i,j,k),
                    f(i,j-1,k) + f(i,j,k),
                    f(i,j+1,k) + f(i,j,k),
                    f(i,j,k-1) + f(i,j,k),
                    f(i,j,k+1) + f(i,j,k)
                );
        }
    }

    return tInterp;
}

}

}

}
