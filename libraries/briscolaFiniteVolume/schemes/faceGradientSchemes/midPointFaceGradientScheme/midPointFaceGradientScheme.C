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

    Grad = Zero;

    const meshField<faceScalar,MeshType>& fd =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    forAllDirections(Grad, d, i, j, k)
        Grad(d,i,j,k) =
            FaceSpace<Type>
            (
                field(d,i-1,j,  k  ) - field(d,i,j,k),
                field(d,i+1,j,  k  ) - field(d,i,j,k),
                field(d,i,  j-1,k  ) - field(d,i,j,k),
                field(d,i,  j+1,k  ) - field(d,i,j,k),
                field(d,i,  j,  k-1) - field(d,i,j,k),
                field(d,i,  j,  k+1) - field(d,i,j,k)
            )
          * fd(d,i,j,k);

    return tGrad;
}

}

}

}
