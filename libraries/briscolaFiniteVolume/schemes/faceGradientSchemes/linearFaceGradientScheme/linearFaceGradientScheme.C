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
tmp<meshField<LowerFaceSpace<Type>,MeshType>>
linearFaceGradientScheme<Type,MeshType>::faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<LowerFaceSpace<Type>,MeshType>> tGrad
    (
        new meshField<LowerFaceSpace<Type>,MeshType>
        (
            "faceGrad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<LowerFaceSpace<Type>,MeshType>& Grad = tGrad.ref();

    Grad = Zero;

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    // The face gradient is defined along the outward normal, for consistency
    // with the flux.

    forAllFaces(Grad, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Grad(d,ijk)[fd] = (field(d,nei) - field(d,ijk)) * delta(d,ijk)[fd*2];
    }

    return tGrad;
}

}

}

}
