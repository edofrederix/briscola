#include "midPointFaceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointFaceGradientScheme<Type,MeshType>::midPointFaceGradientScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointFaceGradientScheme<Type,MeshType>::midPointFaceGradientScheme
(
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
midPointFaceGradientScheme<Type,MeshType>::faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tGrad
    (
        new meshField<FaceSpace<Type>,MeshType>
        (
            "faceGrad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,MeshType>& Grad = tGrad.ref();

    forAll(field, l)
    forAll(field[l], d)
    {
        const meshDirection<faceScalar,MeshType>& fd =
            field.fvMsh().template
            metrics<MeshType>().faceDeltas()[l][d];

        meshDirection<FaceSpace<Type>,MeshType>& G = Grad[l][d];
        const meshDirection<Type,MeshType>& f = field[l][d];

        G.initGhosts();

        forAllCells(f, i, j, k)
        {
            G(i,j,k) =
                FaceSpace<Type>
                (
                    f(i,j,k)-f(i-1,j,  k  ),
                  - f(i,j,k)+f(i+1,j,  k  ),
                    f(i,j,k)-f(i,  j-1,k  ),
                  - f(i,j,k)+f(i,  j+1,k  ),
                    f(i,j,k)-f(i,  j,  k-1),
                  - f(i,j,k)+f(i,  j,  k+1)
                )
              * fd(i,j,k);
        }
    }

    return tGrad;
}

}

}

}
