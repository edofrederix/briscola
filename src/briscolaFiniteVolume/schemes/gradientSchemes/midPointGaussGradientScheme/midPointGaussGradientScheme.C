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
    typedef typename gradientScheme<Type,MeshType>::GradType GType;

    tmp<meshField<GType,MeshType>> tGrad =
        meshField<GType,MeshType>::New
        (
            "grad("+field.name()+")",
            field.fvMsh()
        );

    meshField<GType,MeshType>& Grad = tGrad.ref();

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

        const GType value =
            0.5*fan[fd](d,ijk)*(field(d,ijk) + field(d,nei));

        Grad(d,ijk) += value;
        Grad(d,nei) -= value;
    }

    Grad *= icv;

    return tGrad;
}

}

}

}
