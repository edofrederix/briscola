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
template<template<class> class OpType>
void stencilLinearGaussLaplacianScheme<Type,MeshType>::imLaplacian
(
    linearSystem<stencil,Type,MeshType>& sys,
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .restrict();

    meshField<stencil,MeshType>& A = sys.A();
    meshField<Type,MeshType>& b = sys.b();

    const meshField<faceScalar,MeshType>& fa =
        field.fvMsh().template metrics<MeshType>().faceAreas();

    const meshField<faceScalar,MeshType>& delta =
        field.fvMsh().template metrics<MeshType>().faceDeltas();

    if (isEqOp<OpType>() || isEqMinusOp<OpType>())
    {
        A = Zero;
        b = Zero;
    }

    const OpType<scalar> bop;

    forAllFaces(A, l, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        scalar value = factor*fa(l,d,ijk)[fd*2]*delta(l,d,ijk)[fd*2];

        if (lambdaPtr)
            value *= lambdaPtr->operator()(l,d,ijk)[fd];

        bop(A(l,d,ijk)[fd*2+1], value);
        bop(A(l,d,nei)[fd*2+2], value);

        if (isEqOp<OpType>())
        {
            A(l,d,ijk)[0] -= value;
            A(l,d,nei)[0] -= value;
        }
        else if (isEqMinusOp<OpType>())
        {
            A(l,d,ijk)[0] += value;
            A(l,d,nei)[0] += value;
        }
        else
        {
            bop(A(l,d,ijk)[0], -value);
            bop(A(l,d,nei)[0], -value);
        }
    }

    if (lambdaPtr)
        const_cast<meshField<lowerFaceScalar,MeshType>&>(*lambdaPtr)
       .makeShallow();
}

}

}

}
