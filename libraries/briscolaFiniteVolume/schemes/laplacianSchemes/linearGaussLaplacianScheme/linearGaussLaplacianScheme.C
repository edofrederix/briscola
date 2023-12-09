#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
linearGaussLaplacianScheme<Type,MeshType>::linearGaussLaplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    laplacianScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
linearGaussLaplacianScheme<Type,MeshType>::linearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    laplacianScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    meshField<Type,MeshType>& field
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr).restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>(field)
    );

    linearSystem<stencil,Type,MeshType>& Sys = tSys.ref();

    meshField<stencil,MeshType>& A = Sys.A();

    A = Zero;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = fa(l,d,ijk)[fd*2]*delta(l,d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(l,d,ijk)[fd];

        A(l,d,ijk)[fd*2+1] = value;
        A(l,d,nei)[fd*2+2] = value;

        A(l,d,ijk).center() -= value;
        A(l,d,nei).center() -= value;
    }

    Sys.b() = Zero;

    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussLaplacianScheme<Type,MeshType>::exLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tLap
    (
        new meshField<Type,MeshType>
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,MeshType>& Lap = tLap.ref();

    Lap = Zero;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllFaces(Lap, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = fa(d,ijk)[fd*2]*delta(d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(d,ijk)[fd];

        Lap(d,ijk) += value*(field(d,nei) - field(d,ijk))/cv(d,ijk);
        Lap(d,nei) += value*(field(d,ijk) - field(d,nei))/cv(d,nei);
    }

    return tLap;
}

}

}

}
