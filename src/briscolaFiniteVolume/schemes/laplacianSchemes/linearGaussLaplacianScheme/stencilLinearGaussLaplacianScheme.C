#include "stencilLinearGaussLaplacianScheme.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
stencilLinearGaussLaplacianScheme<Type,MeshType>::
stencilLinearGaussLaplacianScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    linearGaussLaplacianScheme<stencil,Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
stencilLinearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    const faceField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    if (lambdaPtr)
        restrict(*lambdaPtr);

    tmp<linearSystem<stencil,Type,MeshType>> tSys =
        linearSystem<stencil,Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            const_cast<meshField<Type,MeshType>&>(field)
        );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();

    const faceField<scalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const faceField<scalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    A = Zero;

    forAllFaces(delta, fd, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        scalar value = fa[fd](l,d,ijk)*delta[fd](l,d,ijk);

        if (lambdaPtr)
            value *= lambdaPtr->operator[](fd)(l,d,ijk);

        A(l,d,ijk)[fd*2+1] = value;
        A(l,d,nei)[fd*2+2] = value;

        A(l,d,ijk)[0] -= value;
        A(l,d,nei)[0] -= value;
    }

    sys.b() = Zero;

    if (lambdaPtr)
        collapse(*lambdaPtr);

    return tSys;
}

}

}

}
