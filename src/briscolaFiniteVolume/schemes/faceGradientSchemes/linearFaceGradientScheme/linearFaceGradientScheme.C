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
tmp<meshField<FaceSpace<Type>,MeshType>>
linearFaceGradientScheme<Type,MeshType>::faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tGrad =
        meshField<FaceSpace<Type>,MeshType>::New
        (
            "faceGrad("+field.name()+")",
            field.fvMsh()
        );

    meshField<FaceSpace<Type>,MeshType>& Grad = tGrad.ref();

    const faceField<scalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().soa().faceDeltas();

    // The face gradient is defined along the outward normal, for consistency
    // with the flux.

    forAllFaces(Grad, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Grad(d,ijk)[fd*2  ] = (field(d,nei) - field(d,ijk))*delta[fd](d,ijk);
        Grad(d,nei)[fd*2+1] = -Grad(d,ijk)[fd*2];
    }

    return tGrad;
}

}

}

}
