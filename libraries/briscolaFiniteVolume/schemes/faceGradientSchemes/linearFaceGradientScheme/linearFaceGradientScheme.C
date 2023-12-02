#include "linearFaceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearFaceGradientScheme<Type,MeshType>::linearFaceGradientScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearFaceGradientScheme<Type,MeshType>::linearFaceGradientScheme
(
    const fvMesh& fvMsh
)
:
    faceGradientScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
linearFaceGradientScheme<Type,MeshType>::faceGrad
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
        for (int f = 0; f < 6; f++)
            Grad(d,i,j,k)[f] =
                (
                    field(d,labelVector(i,j,k)+faceOffsets[f])
                  - field(d,i,j,k)
                )
              * fd(d,i,j,k)[f];

    return tGrad;
}

}

}

}
