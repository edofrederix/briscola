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
    const fvMesh& fvMsh,
    Istream& is
)
:
    gradientScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<meshField<typename gradientScheme<Type,MeshType>::GradType, MeshType>>
linearGaussGradientScheme<Type,MeshType>::grad
(
    const meshField<Type,MeshType>& field
)
{
    typedef typename gradientScheme<Type,MeshType>::GradType GradType;

    tmp<meshField<GradType,MeshType>> tGrad =
        meshField<GradType,MeshType>::New
        (
            "grad("+field.name()+")",
            field.fvMsh()
        );

    meshField<GradType,MeshType>& Grad = tGrad.ref();

    const meshField<faceVector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();


    forAllCells(Grad, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Grad(d,i,j,k) = Zero;
        #endif

        for (int f = 0; f < 6; f++)
            Grad(d,i,j,k) +=
                fan(d,i,j,k)[f]
              * (
                    fwc(d,i,j,k)[f]*field(d,i,j,k)
                  + fwn(d,i,j,k)[f]*field(d,neighbor(i,j,k,f))
                );

        Grad(d,i,j,k) *= icv(d,i,j,k);
    }

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<Type,staggered>>
linearGaussGradientScheme<Type,MeshType>::stagGrad
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<Type,staggered>> tGrad =
        meshField<Type,staggered>::New
        (
            "stagGrad("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,staggered>& Grad = tGrad.ref();

    const meshField<faceScalar,colocated>& fd =
            field.fvMsh().template metrics<colocated>().faceDeltas();

    forAllCells(Grad, d, i, j, k)
        Grad(d,i,j,k) =
            (
                field(i,j,k)
              - field(lowerNeighbor(i,j,k,d))
            )
          * fd(i,j,k)[d*2];

    return tGrad;
}

}

}

}
