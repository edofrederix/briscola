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

    forAllLevels(Grad, l, d, i, j, k)
        Grad(l,d,i,j,k) =
            0.5
          * (
                (field(l,d,i,j,k)+field(l,d,i-1,j,k))*fan(l,d,i,j,k).left()
              + (field(l,d,i,j,k)+field(l,d,i+1,j,k))*fan(l,d,i,j,k).right()
              + (field(l,d,i,j,k)+field(l,d,i,j-1,k))*fan(l,d,i,j,k).bottom()
              + (field(l,d,i,j,k)+field(l,d,i,j+1,k))*fan(l,d,i,j,k).top()
              + (field(l,d,i,j,k)+field(l,d,i,j,k-1))*fan(l,d,i,j,k).aft()
              + (field(l,d,i,j,k)+field(l,d,i,j,k+1))*fan(l,d,i,j,k).fore()
            )
          / cv(l,d,i,j,k);

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

    forAll(Grad, l)
    {
        forAll(Grad[l], d)
        {
            const labelVector o(staggered::padding[d]);
            const label fo(faceNumber(o));

            forAllCells(Grad[l][d], i, j, k)
            {
                const labelVector ijk(i,j,k);

                Grad(l,d,ijk) =
                    (field(l,0,ijk)-field(l,0,ijk-o))*fd(l,0,ijk)[fo];
            }
        }
    }

    return tGrad;
}

}

}

}
