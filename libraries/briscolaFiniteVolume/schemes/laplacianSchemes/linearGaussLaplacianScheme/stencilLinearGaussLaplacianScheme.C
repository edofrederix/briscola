#include "stencilLinearGaussLaplacianScheme.H"

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
    const meshField<faceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    if (lambdaPtr)
        lambdaPtr->restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>
        (
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    forAllCells(A, l, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        A(l,d,i,j,k)[0] = Zero;
        #endif

        for (int f = 0; f < 6; f++)
        {
            scalar coeff =
                factor*fa(l,d,i,j,k)[f]*delta(l,d,i,j,k)[f];

            if (lambdaPtr)
                coeff *= lambdaPtr->operator()(l,d,i,j,k)[f];

            A(l,d,i,j,k)[f+1] = coeff;
            A(l,d,i,j,k)[0] -= coeff;
        }
    }

    #ifdef NO_BLOCK_ZERO_INIT
    meshField<Type,MeshType>& b = sys.b();
    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) = Zero;
    #endif

    if (lambdaPtr)
        lambdaPtr->makeShallow();

    return tSys;
}

}

}

}
