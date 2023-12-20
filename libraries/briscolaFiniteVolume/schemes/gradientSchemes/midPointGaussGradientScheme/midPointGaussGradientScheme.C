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
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    gradientScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
midPointGaussGradientScheme<Type,MeshType>::midPointGaussGradientScheme
(
    const fvMesh& fvMsh
)
:
    gradientScheme<Type,MeshType>(dictionary(),fvMsh)
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
        Grad(d,i,j,k) =
            0.5
          * (
                (field(d,i,j,k) + field(d,i-1,j,k)) * fan(d,i,j,k).left()
              + (field(d,i,j,k) + field(d,i+1,j,k)) * fan(d,i,j,k).right()
              + (field(d,i,j,k) + field(d,i,j-1,k)) * fan(d,i,j,k).bottom()
              + (field(d,i,j,k) + field(d,i,j+1,k)) * fan(d,i,j,k).top()
              + (field(d,i,j,k) + field(d,i,j,k-1)) * fan(d,i,j,k).aft()
              + (field(d,i,j,k) + field(d,i,j,k+1)) * fan(d,i,j,k).fore()
            )
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
    {
        const labelVector& o = staggered::padding[d];
        const labelVector ijk(i,j,k);

        Grad(d,ijk) =
            (field(0,ijk)-field(0,ijk-o))*fd(0,ijk)[d*2];
    }

    return tGrad;
}

}

}

}
