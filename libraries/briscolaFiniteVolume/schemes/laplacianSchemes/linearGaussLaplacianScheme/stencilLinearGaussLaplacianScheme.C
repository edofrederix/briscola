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
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<stencil,Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
stencilLinearGaussLaplacianScheme<Type,MeshType>::
stencilLinearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<stencil,Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
stencilLinearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>
        (
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    A = Zero;
    b = Zero;

    forAllCells(A, l, d, i, j, k)
    {
        labelVector ijk(i,j,k);

        for (int fd = 0; fd < 3; fd++)
        {
            labelVector nei(upperNei(ijk,fd));

            scalar lower = factor*fa(l,d,ijk)[fd*2  ]*delta(l,d,ijk)[fd*2  ];
            scalar upper = factor*fa(l,d,ijk)[fd*2+1]*delta(l,d,ijk)[fd*2+1];

            if (lambdaPtr)
            {
                lower *= lambdaPtr->operator()(l,d,ijk)[fd];
                upper *= lambdaPtr->operator()(l,d,nei)[fd];
            }

            A(l,d,ijk)[2*fd+1] = lower;
            A(l,d,ijk)[2*fd+2] = upper;

            A(l,d,ijk)[0] -= lower + upper;
        }
    }

    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .makeShallow();

    return tSys;
}

}

}

}
