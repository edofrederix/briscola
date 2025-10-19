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

    tmp<meshField<GradType,MeshType>> tGrad =
        meshField<GradType,MeshType>::New
        (
            "grad("+field.name()+")",
            field.fvMsh()
        );

    meshField<GradType,MeshType>& Grad = tGrad.ref();

    const faceField<vector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    #ifdef NO_BLOCK_ZERO_INIT
    Grad = Zero;
    #endif

    forAllFaces(fan, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const GradType value =
            0.5*fan[fd](d,ijk)*(field(d,ijk) + field(d,nei));

        Grad(d,ijk) += value;
        Grad(d,nei) -= value;
    }

    Grad *= icv;

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<Type,staggered>>
midPointGaussGradientScheme<Type,MeshType>::stagGrad
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

    meshField<Type,staggered>& grad = tGrad.ref();

    const faceField<scalar,colocated>& delta =
            field.fvMsh().template metrics<colocated>().faceDeltas();

    forAllCells(grad, d, i, j, k)
        grad(d,i,j,k) =
            (field(i,j,k) - field(lowerNeighbor(i,j,k,d)))*delta[d](i,j,k);

    return tGrad;
}

}

}

}
