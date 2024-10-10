#include "midPointGaussGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
midPointGaussGradientScheme<Type,MeshType>::midPointGaussGradientScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    gradientScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<meshField<typename gradientScheme<Type,MeshType>::GradType, MeshType>>
midPointGaussGradientScheme<Type,MeshType>::grad
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

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllCells(Grad, d, i, j, k)
        for (int f = 0; f < 6; f++)
            Grad(d,i,j,k) +=
                0.5
              * (
                    field(d,i,j,k)
                  + field(d,nei(i,j,k,f))
                )
              * fan(d,i,j,k)[f]
              / cv(d,i,j,k);

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<Type,staggered>>
midPointGaussGradientScheme<Type,MeshType>::stagGrad
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
              - field(lowerNei(i,j,k,d))
            )
          * fd(i,j,k)[d*2];

    return tGrad;
}

}

}

}
