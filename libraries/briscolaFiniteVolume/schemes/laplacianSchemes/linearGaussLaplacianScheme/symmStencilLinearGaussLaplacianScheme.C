#include "symmStencilLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
symmStencilLinearGaussLaplacianScheme<Type,MeshType>::
symmStencilLinearGaussLaplacianScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    linearGaussLaplacianScheme<symmStencil,Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<symmStencil,Type,MeshType>>
symmStencilLinearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .restrict();

    tmp<linearSystem<symmStencil,Type,MeshType>> tSys
    (
        new linearSystem<symmStencil,Type,MeshType>
        (
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<symmStencil,Type,MeshType>& sys = tSys.ref();

    meshField<symmStencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    A = Zero;
    b = Zero;

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

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

            A(l,d,ijk)[fd+1] =  lower;
            A(l,d,nei)[fd+1] =  upper;

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
