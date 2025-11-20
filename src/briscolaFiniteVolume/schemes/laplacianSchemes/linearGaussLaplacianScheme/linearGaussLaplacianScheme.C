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
    const fvMesh& fvMsh,
    Istream& is
)
:
    laplacianScheme<stencil,Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
linearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    const faceField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    bool shallow = false;

    if (lambdaPtr && lambdaPtr->shallow())
    {
        shallow = true;
        restrict(*lambdaPtr);
    }

    tmp<linearSystem<stencil,Type,MeshType>> tSys =
        linearSystem<stencil,Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            const_cast<meshField<Type,MeshType>&>(field)
        );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    sys.diagonal() = false;

    meshField<stencil,MeshType>& A = sys.A();

    const faceField<scalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const faceField<scalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    A = Zero;

    if (lambdaPtr)
    {
        forAllFaces(delta, fd, l, d, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(i,j,k,fd));

            scalar value =
                fa[fd](l,d,ijk)*delta[fd](l,d,ijk)
              * lambdaPtr->operator[](fd)(l,d,ijk);

            A(l,d,ijk)[fd*2+1] = value;
            A(l,d,nei)[fd*2+2] = value;

            A(l,d,ijk)[0] -= value;
            A(l,d,nei)[0] -= value;
        }
    }
    else
    {
        forAllFaces(delta, fd, l, d, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(i,j,k,fd));

            scalar value = fa[fd](l,d,ijk)*delta[fd](l,d,ijk);

            A(l,d,ijk)[fd*2+1] = value;
            A(l,d,nei)[fd*2+2] = value;

            A(l,d,ijk)[0] -= value;
            A(l,d,nei)[0] -= value;
        }
    }

    sys.b() = Zero;

    if (lambdaPtr && shallow)
        collapse(*lambdaPtr);

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
linearGaussLaplacianScheme<Type,MeshType>::exLaplacian
(
    const faceField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tLap =
        meshField<Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,MeshType>& Lap = tLap.ref();

    const faceField<scalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const faceField<scalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    Lap = Zero;

    forAllFaces(fa, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value =
            fa[fd](d,ijk)*delta[fd](d,ijk)
          * (field(d,nei) - field(d,ijk))
          * (lambdaPtr ? lambdaPtr->operator[](fd)(d,ijk) : 1.0);

        Lap(d,ijk) += value;
        Lap(d,nei) -= value;
    }

    Lap *= icv;

    return tLap;
}

}

}

}
