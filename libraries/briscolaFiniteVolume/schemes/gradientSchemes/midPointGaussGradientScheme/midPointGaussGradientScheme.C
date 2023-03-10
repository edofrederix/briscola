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

    forAll(field, l)
    forAll(field[l], d)
    {
        const meshDirection<faceVector,MeshType>& fan =
            field.fvMsh().template
            metrics<MeshType>().faceAreaNormals()[l][d];

        const meshDirection<scalar,MeshType>& cv =
            field.fvMsh().template
            metrics<MeshType>().cellVolumes()[l][d];

        meshDirection<GradType,MeshType>& G = Grad[l][d];
        const meshDirection<Type,MeshType>& f = field[l][d];

        G = Zero;

        forAllCells(f, i, j, k)
        {
            G(i,j,k) =
                0.5
              * (
                    (f(i,j,k) + f(i-1,j,k)) * fan(i,j,k).left()
                  + (f(i,j,k) + f(i+1,j,k)) * fan(i,j,k).right()
                  + (f(i,j,k) + f(i,j-1,k)) * fan(i,j,k).bottom()
                  + (f(i,j,k) + f(i,j+1,k)) * fan(i,j,k).top()
                  + (f(i,j,k) + f(i,j,k-1)) * fan(i,j,k).aft()
                  + (f(i,j,k) + f(i,j,k+1)) * fan(i,j,k).fore()
                )
              / cv(i,j,k);
        }
    }

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

    forAll(Grad, l)
    {
        const meshDirection<Type,colocated>& f = field[l][0];

        const meshDirection<faceScalar,colocated>& fd =
            field.fvMsh().template
            metrics<colocated>().faceDeltas()[l][0];

        forAll(Grad[l], d)
        {
            const labelVector o(staggered::padding[d]);
            const label fo(faceNumber(o));

            meshDirection<Type,staggered>& G = Grad[l][d];

            G = Zero;

            forAllCells(f, i, j, k)
            {
                const labelVector ijk(i,j,k);

                G(ijk) = (f(ijk)-f(ijk-o))*fd(ijk)[fo];
            }
        }
    }

    return tGrad;
}

}

}

}
