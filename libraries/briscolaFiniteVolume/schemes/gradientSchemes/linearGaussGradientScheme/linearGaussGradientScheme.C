#include "linearGaussGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussGradientScheme<Type,MeshType>::linearGaussGradientScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    gradientScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearGaussGradientScheme<Type,MeshType>::linearGaussGradientScheme
(
    const fvMesh& fvMsh
)
:
    gradientScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<typename gradientScheme<Type,MeshType>::GradType, MeshType>>
linearGaussGradientScheme<Type,MeshType>::grad
(
    const meshField<Type,MeshType>& field
)
{
    typedef typename gradientScheme<Type,MeshType>::GradType GradType;

    tmp<meshField<GradType,MeshType>> tGrad
    (
        new meshField<GradType,MeshType>
        (
            "grad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<GradType,MeshType>& Grad = tGrad.ref();

    Grad = Zero;

    const meshField<faceVector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllCells(Grad, d, i, j, k)
        for (int f = 0; f < 6; f++)
            Grad(d,i,j,k) +=
                (
                    fwc(d,i,j,k)[f]*field(d,i,j,k)
                  + fwn(d,i,j,k)[f]*field(d,labelVector(i,j,k)+faceOffsets[f])
                )
              * fan(d,i,j,k)[f]
              / cv(d,i,j,k);

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<Type,staggered>>
linearGaussGradientScheme<Type,MeshType>::stagGrad
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<Type,staggered>> tGrad
    (
        new meshField<Type,staggered>
        (
            "stagGrad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,staggered>& Grad = tGrad.ref();

    Grad = Zero;

    const meshField<faceScalar,colocated>& fd =
            field.fvMsh().template metrics<colocated>().faceDeltas();

    forAllCells(Grad, d, i, j, k)
        Grad(d,i,j,k) =
            (
                field(i,j,k)
              - field(labelVector(i,j,k)+faceOffsets[d*2])
            )
          * fd(i,j,k)[d*2];

    return tGrad;
}

}

}

}
