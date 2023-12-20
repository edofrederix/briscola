#include "stagLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
stagLinearGaussLaplacianScheme<Type>::stagLinearGaussLaplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<stencil,Type,staggered>(dict,fvMsh)
{}

template<class Type>
stagLinearGaussLaplacianScheme<Type>::stagLinearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<stencil,Type,staggered>(dictionary(),fvMsh)
{}

template<class Type>
tmp<linearSystem<stencil,Type,staggered>>
stagLinearGaussLaplacianScheme<Type>::imLaplacian
(
    const meshField<lowerFaceScalar,staggered>* lambdaPtr,
    meshField<Type,staggered>& field
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,staggered>&>(*lambdaPtr)
       .restrict();

    tmp<linearSystem<stencil,Type,staggered>> tSys
    (
        new linearSystem<stencil,Type,staggered>(field)
    );

    linearSystem<stencil,Type,staggered>& Sys = tSys.ref();

    meshField<stencil,staggered>& A = Sys.A();

    A = Zero;

    const meshField<faceScalar,staggered>& fa =
        field.fvMsh().template metrics<staggered>().faceAreas();

    const meshField<faceScalar,staggered>& delta =
        field.fvMsh().template metrics<staggered>().faceDeltas();

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = fa(l,d,ijk)[fd*2]*delta(l,d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(l,d,ijk)[fd];

        A(l,d,ijk)[fd*2+1] = value;
        A(l,d,nei)[fd*2+2] = value;

        A(l,d,ijk)[0] -= value;
        A(l,d,nei)[0] -= value;
    }

    Sys.b() = Zero;

    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,staggered>&>(*lambdaPtr)
       .makeShallow();

    return tSys;
}

}

}

}
