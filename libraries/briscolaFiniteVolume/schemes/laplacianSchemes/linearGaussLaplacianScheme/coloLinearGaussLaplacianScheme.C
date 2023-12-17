#include "coloLinearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
coloLinearGaussLaplacianScheme<Type>::coloLinearGaussLaplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<symmStencil,Type,colocated>(dict,fvMsh)
{}

template<class Type>
coloLinearGaussLaplacianScheme<Type>::coloLinearGaussLaplacianScheme
(
    const fvMesh& fvMsh
)
:
    linearGaussLaplacianScheme<symmStencil,Type,colocated>(dictionary(),fvMsh)
{}

template<class Type>
tmp<linearSystem<symmStencil,Type,colocated>>
coloLinearGaussLaplacianScheme<Type>::imLaplacian
(
    const meshField<lowerFaceScalar,colocated>* lambdaPtr,
    meshField<Type,colocated>& field
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,colocated>&>(*lambdaPtr)
       .restrict();

    tmp<linearSystem<symmStencil,Type,colocated>> tSys
    (
        new linearSystem<symmStencil,Type,colocated>(field)
    );

    linearSystem<symmStencil,Type,colocated>& Sys = tSys.ref();

    meshField<symmStencil,colocated>& A = Sys.A();

    A = Zero;

    const meshField<faceScalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    const meshField<faceScalar,colocated>& delta =
        field.fvMsh().template metrics<colocated>().faceDeltas();

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = fa(l,d,ijk)[fd*2]*delta(l,d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(l,d,ijk)[fd];

        A(l,d,ijk)[fd+1] = value;

        A(l,d,ijk)[0] -= value;
        A(l,d,nei)[0] -= value;
    }

    Sys.b() = Zero;

    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,colocated>&>(*lambdaPtr)
       .makeShallow();

    return tSys;
}

}

}

}
