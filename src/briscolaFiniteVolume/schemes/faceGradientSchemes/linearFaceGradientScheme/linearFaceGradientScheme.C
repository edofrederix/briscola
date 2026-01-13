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
    const fvMesh& fvMsh,
    Istream& is
)
:
    faceGradientScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<faceField<Type,MeshType>>
linearFaceGradientScheme<Type,MeshType>::faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    tmp<faceField<Type,MeshType>> tGrad =
        faceField<Type,MeshType>::New
        (
            "faceGrad("+field.name()+")",
            field.fvMsh()
        );

    faceField<Type,MeshType>& grad = tGrad.ref();

    const faceField<scalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    // The face gradient is defined along the outward normal, for consistency
    // with the flux.

    forAllFaces(grad, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        grad[fd](d,ijk) =
            (field(d,nei) - field(d,ijk))*delta[fd](d,ijk);
    }

    return tGrad;
}

}

}

}
